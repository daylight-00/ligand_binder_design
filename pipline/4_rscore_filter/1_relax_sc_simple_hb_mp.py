import os
import sys
import pandas as pd
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from xml_relax_after_ligMPNN_PTM import XML_BSITE_REPACK_MIN_BETA

def repack_pose(work_pose, nres, packable_res, xml_in, tag, cst_fn):
    import pyrosetta.distributed.tasks.rosetta_scripts as distributed_rosetta_scripts
    repack_protocol = xml_in.format(packable_res, cst_fn)
    return distributed_rosetta_scripts.SingleoutputRosettaScriptsTask(repack_protocol), work_pose

class RelaxScore:
    def __init__(self, pdbfn, target_hb_atms, repack_res=[], param_fn=None):
        self.pdbfn = pdbfn.strip()
        self.target_hb_atms = target_hb_atms.split(',')
        self.param_fn = param_fn
        self.tag = self.pdbfn.split('/')[-1].split('.pdb')[0]
        self.output_dir = '/'.join(self.pdbfn.split('/')[:-1])
        self.cst_fn = f'{self.output_dir}/{self.tag}.cst'
        self.repack_res = repack_res
        self.repack_all = len(self.repack_res) == 0

    def bsite_repack_min(self, pose_work):
        self.nres = pose_work.total_residue()
        if self.repack_all:
            self.repack_res = list(range(1, self.nres))
        repackable_res_str = ','.join(str(resno) for resno in self.repack_res)
        xml_obj, pose_work = repack_pose(pose_work, self.nres, repackable_res_str, XML_BSITE_REPACK_MIN_BETA, self.tag, self.cst_fn)
        xml_obj.setup()
        pose_work = xml_obj(pose_work)
        out_pdb = f'{self.output_dir}/{self.tag}_betarelax.pdb'
        return pose_work, out_pdb

    def calc_hb(self, pose_work):
        import pyrosetta
        import pyrosetta.distributed.packed_pose as packed_pose
        full_pose = packed_pose.to_pose(pose_work)
        hbond_set = pyrosetta.rosetta.core.scoring.hbonds.HBondSet()
        full_pose.update_residue_neighbors()
        pyrosetta.rosetta.core.scoring.hbonds.fill_hbond_set(full_pose, False, hbond_set)
        lig_res = full_pose.residue(self.nres)
        lig_atm_hb = defaultdict(list)

        for lig_atmName in self.target_hb_atms:
            atm_idx = lig_res.atom_index(lig_atmName)
            atm_id = pyrosetta.rosetta.core.id.AtomID(atm_idx, self.nres)
            found_hbs = hbond_set.atom_hbonds(atm_id)
            for hb in found_hbs:
                don_resNo = hb.don_res()
                don_atmName = full_pose.residue(don_resNo).atom_name(hb.don_hatm())
                acc_resNo = hb.acc_res()
                acc_atmName = full_pose.residue(acc_resNo).atom_name(hb.acc_atm())
                hb_atm = {don_resNo: don_atmName, acc_resNo: acc_atmName}
                hb_res_pair = [don_resNo, acc_resNo]
                for i_res, resno in enumerate(hb_res_pair):
                    if resno == self.nres:
                        other_resno = hb_res_pair[1 - i_res]
                        if other_resno == resno:
                            continue
                        lig_atm_hb[hb_atm[resno].strip()].append((other_resno, hb_atm[other_resno]))

        hb_sc = {}
        for lig_atmName in self.target_hb_atms:
            tmp = list(set(lig_atm_hb.get(lig_atmName, [])))
            hb_sc[f'{lig_atmName}_hbond'] = len(tmp)
        return hb_sc

pyrosetta_threads = 2
def process_one_pdb(pdb_path, ligname, lig_parm_fn, target_hbatm_names, dump_pdb, output_csv):
    import pyrosetta
    import pyrosetta.distributed.packed_pose as packed_pose
    from get_pock_res_by_dist_lig import PocketPDB
    from gen_prot_lig_dist_cst2 import extract_dist_cst_from_pdb_use_allatm

    init_args = [
        '-mute', 'all',
        '-beta',
        '-gen_potential',
        '-in:file:native', pdb_path,
        '-extra_res_fa', lig_parm_fn,
        # f'-multithreading:total_threads {pyrosetta_threads}',
        # f'-multithreading:interaction_graph_threads {pyrosetta_threads}'
    ]
    init_flags = ' '.join(init_args)
    pyrosetta.init(init_flags)
    tmp_pdb = PocketPDB(pdb_path, ligname)
    repack_res = tmp_pdb.find_close_lig_contact_lig()
    pdb_score = RelaxScore(pdb_path, target_hbatm_names, repack_res=repack_res)
    csts = extract_dist_cst_from_pdb_use_allatm(pdb_path, ligname)
    with open(pdb_score.cst_fn, 'w') as fout:
        fout.write('\n'.join(csts) + '\n')

    pose_init = pyrosetta.pose_from_pdb(pdb_path)
    pose_work = pose_init.clone()
    pose_work, out_pdb = pdb_score.bsite_repack_min(pose_work)
    sc = pose_work.scores
    sc.update(pdb_score.calc_hb(pose_work))
    if dump_pdb:
        packed_pose.to_pose(pose_work).dump_pdb(out_pdb)
    sc['tag'] = pdb_score.tag

    df = pd.DataFrame.from_records([sc])
    header = not os.path.exists(output_csv)
    df.to_csv(output_csv, mode='a', header=header, index=False)

def main():
    mapping = {
    'ligand_params_path': {
        "lmpnn": {
            "pht": "../0_params/PHT.params",
            "cbz": "..."
        },
        "af3": {
            "pht": "...",
            "cbz": "..."
        },
        "boltz": {
            "pht": "...",
            "cbz": "..."
        },
    },
    'target_atm_for_cst': {
        "lmpnn": {
            "pht": "O1,O2,N1,N2",
            "cbz": "..."
        },
        "af3": {
            "pht": "...",
            "cbz": "..."
        },
        "boltz": {
            "pht": "...",
            "cbz": "..."
        },
    }
}

    pdblist = sys.argv[1]
    # ligname = sys.argv[2]
    # lig_parm_fn = sys.argv[3]
    # target_hbatm_names = sys.argv[4]
    # dump_pdb = len(sys.argv) > 5 and sys.argv[5] == 'dump_pdb'
    dump_pdb = True

    # output_csv = f'betagen.csv'
    # if os.path.exists(output_csv):
    #     os.remove(output_csv)

    df_pdb = pd.read_parquet(pdblist)

    failed = []
    col_1 = 'lmpnn'
    col_2 = 'diffusion'
    with ProcessPoolExecutor() as executor:
        futures = {}
        for idx, row in df_pdb.iterrows():
            pdb_path = row[(col_1, 'path')]
            tag = row[(col_2, 'tag')]
            batch_name = row[(col_2, 'batch')]
            lig_name = batch_name.split('_')[0]
            if 'af3' in tag: model = 'af3'
            elif 'boltz' in tag: model = 'boltz'
            else: model = 'lmpnn'
            lig_parm_fn = mapping['ligand_params_path'][model][lig_name]
            target_hbatm_names = mapping['target_atm_for_cst'][model][lig_name]
            ligname = batch_name.split('_')[0].upper()
            # ligname = 'LIG'
            output_csv = f'betagen_{model}_{lig_name}.csv'
            futures[executor.submit(process_one_pdb, pdb_path, ligname, lig_parm_fn, target_hbatm_names, dump_pdb, output_csv)] = pdb_path
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing PDBs"):
            pdb_path = futures[future]
            try:
                future.result()
            except Exception as e:
                print(f"Error on {pdb_path}: {e}")
                failed.append(pdb_path)

    if failed:
        if os.environ.get('SLURM_ARRAY_TASK_ID'):
            SLURM_ARRAY_TASK_ID = os.environ['SLURM_ARRAY_TASK_ID']
            fail_out = f'failed_list_{SLURM_ARRAY_TASK_ID}.txt'
        else:
            fail_out = 'failed_list.txt'

        with open(fail_out, 'a') as fout:
            fout.write('\n'.join(failed) + '\n')
        print(f"{len(failed)} PDBs failed. Appended to '{fail_out}'.")

if __name__ == '__main__':
    main()

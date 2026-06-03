#!/usr/bin/env python3
import os
import sys
import pandas as pd
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import pyrosetta as py

EXTRA_RES = None

def init_pyrosetta(extra_res_path: str):
    global EXTRA_RES
    EXTRA_RES = extra_res_path
    py.init(f"-mute all -beta -gen_potential -extra_res_fa {extra_res_path}")

def get_lig_atom_sasa(pose) -> dict:
    atm_depth = py.rosetta.core.scoring.atomic_depth.atomic_depth(pose, 3.0)
    tnr = pose.total_residue()
    lig_res = pose.residue(tnr)
    natm = lig_res.natoms()

    result = {}
    for i in range(1, natm + 1):
        atm_name = lig_res.atom_name(i).strip()
        atm_id = py.rosetta.core.id.AtomID(
            lig_res.atom_index(atm_name), tnr
        )
        result[atm_name] = atm_depth(atm_id)

    result["sum"] = sum(result.values())
    return result

def process_pdb(pdbfile):
    base = os.path.basename(pdbfile).rsplit('.', 1)[0]
    try:
        pose = py.pose_from_pdb(pdbfile)
        sasa_dict = get_lig_atom_sasa(pose)
        return {"id": base, **sasa_dict}
    except Exception as e:
        print(f"[Error] {pdbfile}: {e}", file=sys.stderr)
        return None

def main():
    if len(sys.argv) != 3:
        print("Usage: script.py <input_csv> <extra_res_file>")
        sys.exit(1)

    input_csv = sys.argv[1]
    extra_res = sys.argv[2]

    df = pd.read_csv(input_csv)
    pdb_paths = df["path"].tolist()

    n_procs = cpu_count()
    chunksize = max(1, len(pdb_paths) // (n_procs * 4))

    with Pool(
        processes=n_procs,
        initializer=init_pyrosetta,
        initargs=(extra_res,)
    ) as pool:
        iterator = pool.imap_unordered(process_pdb, pdb_paths, chunksize=chunksize)
        results = list(tqdm(
            iterator,
            total=len(pdb_paths),
            desc="Processing PDB files"
        ))

    rows = [r for r in results if r]
    df_res = pd.DataFrame(rows)

    cols = ['id', 'sum'] + [c for c in df_res.columns if c not in ('id', 'sum')]
    df_res = df_res[cols]

    output_csv = "lig_atm_depth.csv"
    df_res.to_csv(output_csv, index=False)
    print(f"Saved results to {output_csv}")

if __name__ == "__main__":
    main()

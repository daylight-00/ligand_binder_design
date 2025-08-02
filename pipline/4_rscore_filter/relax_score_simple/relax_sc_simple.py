import os
import sys
import glob
import numpy as np
from collections import defaultdict
import pandas as pd

import pyrosetta
from pyrosetta.rosetta.protocols.simple_moves import *
import pyrosetta.rosetta.protocols.rosetta_scripts as rosetta_scripts

import pyrosetta.distributed.io as io
import pyrosetta.distributed.packed_pose as packed_pose
import pyrosetta.distributed.tasks.rosetta_scripts as rosetta_scripts
import pyrosetta.distributed.tasks.score as score

from xml_relax_after_ligMPNN_PTM import XML_BSITE_REPACK_MIN_BETA
from gen_prot_lig_dist_cst2 import extract_dist_cst_from_pdb_use_allatm,CST_STDERR
from get_pock_res_by_dist import PocketPDB

def repack_pose(work_pose, nres, packable_res, xml_in, tag, cst_fn):
    #cst_fn is global var.
    repack_protocol = xml_in.format(packable_res,cst_fn)
#    repack_protocol = xml_in.format(packable_res,nres)
#    pyrosetta.rosetta.core.pose.setPoseExtraScore(work_pose,"tag",'%s'%tag)
    return rosetta_scripts.SingleoutputRosettaScriptsTask(repack_protocol), work_pose

class RelaxScore:
    def __init__(self,pdbfn,repack_res=[],param_fn=None):
        self.pdbfn = pdbfn.strip()
#        self.rosetta_ptm_patch = rosetta_ptm_patch
#        self.target_hb_atms = target_hb_atms.split(',')
        self.param_fn = param_fn
        self.tag = self.pdbfn.split('/')[-1].split('.pdb')[0]
        self.output_dir = '/'.join(self.pdbfn.split('/')[:-1])
        self.cst_fn = f'{self.output_dir}/{self.tag}.cst'
        self.repack_res = repack_res
        self.repack_all = False
        if (len(self.repack_res) == 0):
            self.repack_all = True
    def bsite_repack_min(self,pose_work):
        self.nres = len(pose_work)
        #Repack all if not any residue is assigned
        if self.repack_all and len(self.repack_res) == 0:
            self.repack_res = []
            for ires in range(1,self.nres):
                self.repack_res.append(ires)
            #
        repackable_res_str = ','.join(['%d'%resno for resno in self.repack_res])
        xml_obj, pose_work = repack_pose(pose_work,self.nres,repackable_res_str,XML_BSITE_REPACK_MIN_BETA,self.tag,self.cst_fn)
        xml_obj.setup()
        pose_work = xml_obj(pose_work)
        out_pdb = '%s/%s_betarelax.pdb'%(self.output_dir,self.tag)
        return pose_work, out_pdb


def main():
    pdblist = sys.argv[1]
    mod_resname = sys.argv[2] #Name of target residue (modified)
#    target_seq = sys.argv[2]
#    target_hbatm_names = sys.argv[3]
    dump_pdb = False
    if len(sys.argv) > 3 and sys.argv[3] == 'dump_pdb':
        dump_pdb = True
    #
    #
    #TODO: to check if this can be recognized by rosetta
#    ptm_resname = target_seq.split("(")[1].split(")")[0]
#    rosetta_ptm_patch = ptm_patch[ptm_resname]
    #
    df_s = []
    with open(pdblist) as fp:
        for line in fp:
            pyrosetta.init(
                "-beta -gen_potential -extra_res_fa /home/hwjang/aipd/250621/2_/tmp/PHT.params -in:file:native %s" % line.strip()
            )
#            pyrosetta.init('-beta -in:file:native %s'%line.strip())
#            pyrosetta.init('-beta -gen_potential -in:file:native %s'%line.strip())
            #get repackable residues (treating chainB as ligand)
            tmp_pdb = PocketPDB(line.strip())
            repack_res = tmp_pdb.find_close_lig_contact_chainB()
            #
            pdb_score = RelaxScore(line.strip(),repack_res=repack_res)
            #
            csts = extract_dist_cst_from_pdb_use_allatm(line.strip(), mod_resname)
            fout = open(pdb_score.cst_fn,'w')
            fout.write('%s\n'%('\n'.join(csts)))
            fout.close()
            pose_init = pyrosetta.pose_from_pdb(line.strip())
            pose_work = pose_init.clone()
            pose_work, out_pdb = pdb_score.bsite_repack_min(pose_work)
            sc = pose_work.scores
            print (sc)
#            hb_d = pdb_score.calc_hb(pose_work)
#            sc.update(hb_d)
            if dump_pdb:
                packed_pose.to_pose(pose_work).dump_pdb(f'{out_pdb}')
            sc['tag'] = pdb_score.tag
            #
            df_s.append(pd.DataFrame.from_records([sc]))
    #
    sc_df = pd.DataFrame()
    if (len(df_s) > 0):
        sc_df = pd.concat(df_s,ignore_index=True)
 #       hbatm_s = target_hbatm_names.split(',')
 #       hb_tags = []
 #       for atmname in hbatm_s:
 #           hb_tags.append('%s_hbond'%atmname)
 #       sc_df['phos_nhb'] = sc_df[hb_tags].sum(axis=1)
#        sc_df['phos_nhb'] = sc_df['O1P_hbond']+sc_df['O2P_hbond']+sc_df['O3P_hbond']+sc_df['OH_hbond']
    #
    sc_df.to_csv('./%s_betagen.csv'%pdblist.split('/')[-1])
    return sc_df

if __name__ == '__main__':
    main()

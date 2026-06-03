import sys, os
import PLACER
import pandas as pd
from tqdm import tqdm

placer = PLACER.PLACER()

def process(fname, ligand_file, predict_ligand, odir, nsamples=50, rerank="prmsd"):
    if ".cif.gz" in os.path.basename(fname):
        label = os.path.basename(fname).replace(".cif.gz", "")
    elif ".cif" in os.path.basename(fname):
        label = os.path.basename(fname).replace(".cif", "")
    elif ".pdb" in os.path.basename(fname):
        label = os.path.basename(fname).replace(".pdb", "")
    outfile_prefix = odir + "/" + label

    placer_input = PLACER.PLACERinput()
    if fname.endswith(".pdb"):
        placer_input.pdb(fname)
    elif fname.endswith(".cif") or fname.endswith(".cif.gz"):
        placer_input.cif(fname)

    placer_input.name(os.path.basename(fname).replace(".cif", "").replace(".cif.gz", "").replace(".pdb", ""))
    ligand_ref = {predict_ligand: ligand_file}
    placer_input.ligand_reference(ligand_ref)
    outputs = placer.run(placer_input, nsamples)
    os.makedirs(odir, exist_ok=True)
    PLACER.protocol.dump_output(output_dict=outputs, filename=outfile_prefix, rerank=rerank)

nsamples = 50
odir = "output"
rerank = "prmsd"
# predict_ligand = "PHT"
projects = ['pht_demo']
mappings = {
    'pht': '../0_params/Conformer3D_COMPOUND_CID_1775.sdf',
}

input_file = sys.argv[1]
col_1 = 'lmpnn'
df = pd.read_parquet(input_file)
for idx, row in tqdm(df.iterrows(), total=len(df)):
    fname = row[(col_1, 'path')]
    batch_name = row[('diffusion', 'batch')]
    ligname = batch_name.split('_')[0]
    ligand_file = mappings[ligname.lower()]
    print(fname, ligand_file, ligname.upper())
    process(fname, ligand_file, ligname.upper(), odir, nsamples, rerank)

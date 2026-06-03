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

def get_dssp_string(pdbfile: str) -> str:
    pose = py.pose_from_pdb(pdbfile)
    mover = py.rosetta.protocols.moves.DsspMover()
    mover.apply(pose)
    return pose.secstruct()

def process_pdb(pdbfile: str):
    base = os.path.basename(pdbfile).rsplit('.', 1)[0]
    try:
        ss = get_dssp_string(pdbfile)
        return {"id": base, "dssp": ss}
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
            desc="Processing DSSP"
        ))

    rows = [r for r in results if r]
    df_res = pd.DataFrame(rows)
    df_res = df_res[['id', 'dssp']]

    output_csv = "dssp.csv"
    df_res.to_csv(output_csv, index=False)
    print(f"Saved DSSP results to {output_csv}")

if __name__ == "__main__":
    main()

# Simple Pipeline for Ligand Binder Design

**Author:** David Hyunyoo Jang  
**Affiliation:** [Artificial Intelligence Protein Design Lab](https://sites.google.com/view/aipdlab)  
**Date:** July 2025

## Overview
This repository provides a simple pipeline for ligand binder design using diffusion models. The pipeline integrates multiple state-of-the-art tools to generate and validate protein binders for small molecule ligands.

## Pipeline Steps

The pipeline consists of 9 main steps:

1. **Parameters Setup** (`0_params/`) - Ligand preparation and parameter files
2. **Backbone Generation** (`1_diffusion/`) - Structure generation using RFDiffusion
3. **Backbone Filtering** (`2_backbone_filter/`) - DSSP and SASA-based filtering
4. **Sequence Design** (`3_lmpnn/`) - Sequence generation using LigandMPNN
5. **Rosetta Scoring** (`4_rscore_filter/`) - Energy-based filtering
6. **AlphaFold3 Prediction** (`5_af3/`) - Structure prediction and validation
7. **Boltz Prediction** (`6_boltz/`) - Alternative structure prediction
8. **PLACER Analysis** (`7_placer/`) - Binding site analysis
9. **RMSD Filtering** (`8_rmsd_filter/`) - Final structure validation

## Features

- **Automated Pipeline**: End-to-end workflow from ligand input to validated binders
- **Multiple Validation Steps**: Combines geometric, energetic, and structural filters
- **Modern AI Tools**: Utilizes RFDiffusion, LigandMPNN, AlphaFold3, and Boltz
- **Scalable**: Designed for SLURM-based cluster environments
- **Performance Optimized**: Multiprocessing enabled for DSSP/SASA/Rosetta scoring and RMSD calculations
- **Enhanced RMSD Analysis**: Biopython-based structure handling with RDKit for ligand symmetry-aware RMSD calculations

## Requirements

- RFDiffusion All-Atom
- LigandMPNN
- Rosetta
- AlphaFold3
- Boltz
- PLACER
- TMalign
- PyMOL
- Python 3.8+
- PyArrow (for parquet file handling)

**Note**: This pipeline uses parquet file format to handle double headers efficiently. Parquet files can be viewed and analyzed using VSCode's [Data Wrangler](https://marketplace.visualstudio.com/items?itemName=ms-toolsai.datawrangler) extension.

## Usage

1. Place your ligand files in the `0_params/` directory
2. Run each pipeline step sequentially using the provided `run.sh` scripts
3. Results will be filtered at each step, with final candidates in `8_rmsd_filter/`

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Gyu Rie Lee - Original concept and supervision

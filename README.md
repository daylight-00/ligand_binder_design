# Reproducible Small-Molecule Binder Design Pipeline with Multi-Model Structural Triage

This repository organizes an end-to-end computational workflow for de novo small-molecule binder design. It connects ligand parameterization, diffusion-based backbone generation, LigandMPNN sequence design, Rosetta scoring, AF3/Boltz-2 refolding, PLACER analysis, and RMSD-based structural triage into a reproducible design pipeline.

The workflow is designed for small-molecule binder discovery campaigns where many generated candidates must be filtered through multiple computational criteria before experimental consideration.

This repository emphasizes:

- reproducible execution of each design stage;
- structured storage of candidate metadata and metrics;
- multi-model structural validation rather than reliance on a single predictor;
- target-aware RMSD calculations for chemically modified ligands;
- scalable execution on SLURM-based computing environments.

## Workflow Overview

The main workflow is organized under `pipeline/` and follows nine sequential stages.

1. **Parameters Setup** (`0_params/`) - Ligand preparation and Rosetta parameter generation
2. **Backbone Sampling** (`1_diffusion/`) - Diffusion-based backbone generation
3. **Backbone Filtering** (`2_backbone_filter/`) - DSSP, SASA, and backbone-quality filtering
4. **Sequence Sampling** (`3_lmpnn/`) - LigandMPNN sequence design
5. **Rosetta Scoring** (`4_rscore_filter/`) - Rosetta scoring and interface-energy filtering
6. **AlphaFold3 Refolding** (`5_af3/`) - AlphaFold3 refolding and structural validation
7. **Boltz Refolding** (`6_boltz/`) - Boltz-2 refolding and alternative model validation
8. **PLACER Analysis** (`7_placer/`) - PLACER-based binding-site analysis
9. **RMSD Filtering** (`8_rmsd_filter/`) - Binder, motif, and ligand RMSD filtering

Each stage can be run independently, but the intended use is sequential: outputs from one stage are converted into inputs and metadata tables for the next stage.

## What This Repository Adds

This repository does not claim to introduce a new protein design model. Its value is in workflow integration and research software engineering for practical binder design campaigns.

### 1. Tutorial-to-pipeline conversion

The initial workflow was learned from Prof. Lee’s ligand binder design tutorial. This repository converts the tutorial-style sequence of commands and scripts into a reusable pipeline with more consistent inputs, outputs, and stage organization.

### 2. Target-specific design campaigns

The workflow was applied to commissioned small-molecule sensor targets, including phenytoin and carbamazepine. These targets motivated additional handling for ligand geometry, RMSD evaluation, and downstream triage.

### 3. Structured metric lineage

Candidate information is tracked in table-oriented formats, including parquet, so that generation metadata and downstream metrics can be joined across stages.

Typical tracked information includes:

- target and ligand identifiers;
- generation run IDs and input structures;
- backbone filtering metrics;
- LigandMPNN sequence-design outputs;
- Rosetta interface and total-score metrics;
- AF3/Boltz-2 confidence and structural metrics;
- PLACER outputs;
- TMalign and RMSD-based validation metrics;

The goal is to move from a file-centric design workflow to a table-centric design database.

### 4. Multi-model structural triage

The workflow supports structural validation using multiple predictors and scoring tools rather than relying on a single model. In the current campaign structure, AF3 and Boltz-2 are used for refolding and validation, while Rosetta, PLACER, TMalign, and RMSD-based filters are used to evaluate interface quality and structural consistency.

### 5. Symmetry-aware and linker-agnostic ligand RMSD

For small-molecule sensor design, ligand RMSD can be misleading if molecular symmetry or linker modifications are ignored. This repository includes RMSD utilities designed to support:

- symmetry-aware ligand RMSD;
- ligand RMSD after binder-side superposition;
- RMSD calculations that can ignore linker atoms when comparing modified and unmodified ligand representations.

This is especially relevant when the designed binder is intended for sensor development and the experimental ligand representation may include a linker.

### 6. Script modernization and scalability

Several lab-provided scripts were refactored or extended for more robust campaign-scale use. Updates include:

- more uniform command-line and file interfaces;
- multiprocessing for PyRosetta-based scoring and filtering;
- updated structure parsing with modern Python libraries where appropriate;
- Biopython/RDKit-based ligand and structure handling for RMSD workflows;
- parquet-based data handling for large candidate sets.

## Repository layout

- **`main`** — pipeline code, ligand-setup inputs, and phenytoin example run (final filtered table, joined metric-lineage table, and example RMSD structures).
- **`archive/example-run`** — the same phenytoin example run, with all of its outputs, prediction results, and run logs.

## Requirements

The full workflow depends on several external modeling and scoring tools. Depending on which stages are used, requirements may include:

- RFdiffusion / RFdiffusion All-Atom
- LigandMPNN
- PyRosetta
- AlphaFold3
- Boltz-2
- PLACER
- TMalign
- PyMOL
- RDKit
- Biopython
- Pandas (with PyArrow for parquet handling)

**Note**: This workflow uses parquet file format to handle double headers efficiently. Parquet files can be viewed and analyzed using VSCode's [Data Wrangler](https://marketplace.visualstudio.com/items?itemName=ms-toolsai.datawrangler) extension.

## Notes on RMSD Evaluation

The RMSD stage distinguishes several structural questions:

- **Global binder RMSD**: Does the designed binder fold similarly after prediction?
- **Local motif RMSD**: Are binding-site or motif residues preserved?
- **Ligand RMSD**: Does the ligand remain in the intended binding pose after binder-side alignment?
- **Symmetry-aware ligand RMSD**: Are chemically equivalent atom mappings handled correctly?
- **Linker-agnostic ligand RMSD**: Can ligand pose be evaluated without penalizing added linker atoms?

This is important because a ligand may be represented differently between design, modeling, and sensor-development contexts.

## Relationship to Other Projects

This repository is complementary to the [*ligandMPNN_FR*](https://github.com/daylight-00/ligandMPNN_FR) project.

```text
ligand_binder_design:
    end-to-end target-design workflow
    generation → filtering → sequence design → refolding → metric lineage

ligandMPNN_FR:
    metric-guided LigandMPNN–FastRelax recycling
    optimization of backbone-conditioned sequence sampling spaces
```

Together, the two projects represent different layers of the small-molecule binder design process: a practical binder-design campaign pipeline and a more focused optimization strategy for improving candidate pools.

## Acknowledgements

This project was developed during a research internship in the [*Artificial Intelligence Protein Design Lab*](https://sites.google.com/view/aipdlab) under the supervision of Prof. Gyu Rie Lee. The workflow was initially motivated by Prof. Lee’s binder design tutorial and was extended during target-specific small-molecule binder design campaigns.

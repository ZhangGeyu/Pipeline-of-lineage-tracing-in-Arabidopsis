# Pipeline-of-lineage-tracing-in-Arabidopsis

This repository contains the analysis pipeline for lineage tracing in Arabidopsis, integrating UMIC-seq and TA-clone sequencing data to identify bona fide somatic mutations, reconstruct lineage relationships, and analyze mutation frequencies across samples.

This version (Analysis-Pipeline-v2.0) includes updated workflows, environment definitions, and documentation for reproducible analysis.

## Repository Structure

├── Analysis-pipelines/           # All analysis scripts and workflows
├── Data/                         # Raw and intermediate data used by the pipeline
├── Envs/                         # Conda environment .yml files
├── Figures-and-related-code/     # Scripts and files needed to generate figures
└── README.md                     # This file

### Computational Environments
Conda environment files required to run the analysis are provided in the Envs/ folder. Use these to reproduce the software environment used in this study.
#### UMIC-seq.yml
Environment for UMIC-seq preprocessing, UMI extraction, and clustering.
#### GATK.yml
Environment for whole-genome sequencing mutation calling.
#### Analysis-Pipeline.yml
Environment for downstream scripts including mutation filtering, frequency calculation, and lineage analyses.

#### To create and activate an environment:
conda env create -f Envs/UMIC-seq.yml

conda activate UMIC-seq

### Overview of Pipeline Steps

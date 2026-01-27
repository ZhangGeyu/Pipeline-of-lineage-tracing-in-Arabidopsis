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
#### 1. UMIC-seq Preprocessing & Clustering
Extract sample barcodes and UMIs.

Demultiplex reads by sample.

Perform UMI clustering to group reads with identical or similar UMIs.

Output clustered FASTA for downstream analysis.

#### 2. UMIC-seq Mutation Calling

Merge clustered reads from all samples.

Align merged reads to the reference genome.

Split alignments by readout.

Call mutations using mpileup files.

Filter bona fide mutations.

#### 3. Progeny TA-clone Analysis

Combine forward and reverse Sanger sequencing into consensus sequences.

Map progeny consensus sequences to reference.

Extract progeny mutations.


#### 4. Frequency Analysis & Lineage Inference

Extract readouts for selected copies.

Compute mutation frequencies in parental and progeny samples.

Generate lineage trees and statistical analysis.







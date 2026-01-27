# Pipeline-of-lineage-tracing-in-Arabidopsis

This repository contains the analysis pipeline for lineage tracing in Arabidopsis, integrating UMIC-seq and TA-clone sequencing data to identify bona fide somatic mutations, reconstruct lineage relationships, and analyze mutation frequencies across samples.

This version (Analysis-Pipeline-v2.0) includes updated workflows, environment definitions, and documentation for reproducible analysis.

## Repository Structure

Pipeline-of-lineage-tracing-in-Arabidopsis/

├── Analysis-pipelines/           # Analysis scripts and workflows for lineage tracing

├── Data/                         # Raw and intermediate data used in the analyses

├── Envs/                         # Conda environment files for reproducible computation

├── Figures-and-related-code/     # Scripts and auxiliary files for figure generation

├── LICENSE                       # MIT License for this repository

└── README.md                     # Project overview and usage instructions

## Computational Environments
Conda environment files required to run the analysis are provided in the Envs/ folder. Use these to reproduce the software environment used in this study.
#### 1. UMIC-seq.yml
Environment for UMIC-seq preprocessing, UMI extraction, and clustering.
#### 2. GATK.yml
Environment for whole-genome sequencing mutation calling.
#### 3. Analysis-Pipeline.yml
Environment for downstream scripts including mutation filtering, frequency calculation, and lineage analyses.

#### To create and activate an environment:
conda env create -f Envs/UMIC-seq.yml

conda activate UMIC-seq

## Overview of analysis-pipelines
Codes required to run the analysis are provided in the Analysis-pipelines/ folder.
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

## Raw and intermediate data

Raw and intermediate data files are provided in the Data/ folder

This includes:

#### Barcodes and probe sequences used for UMIC-seq preprocessing.
#### Reference sequences for each plant.
#### Sample metadata files linking barcodes to sample IDs.
#### Mutation information tables for parental and progeny analyses.

## Figures and Visualization
All scripts and auxiliary files used to generate figures in the associated manuscript are in Figures-and-related-code/

## Citation
If you use this pipeline in your research, please cite the corresponding manuscript and include the version tag (e.g., v2.0) from this repository.

## License
This project is licensed under the MIT License.

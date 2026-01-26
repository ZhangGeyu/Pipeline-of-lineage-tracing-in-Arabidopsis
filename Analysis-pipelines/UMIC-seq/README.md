This directory contains the data processing pipeline for UMIC-seq experiments.
The pipeline is designed to identify bona fide somatic mutations by tracking UMI-tagged reads across Arabidopsis samples.

The workflow below uses the Plant1 Leaf UMIC-seq library as an example.
The corresponding input files can be found in: Pipeline-of-lineage-tracing-in-Arabidopsis/Data/UMIC-seq

# Pipeline Steps
## 1. Extraction of Sample Barcodes (Sample_BC)

In this step, sample barcodes are extracted from raw sequencing reads.
These barcodes are used to assign and group reads according to their originating samples.

### Step1.1: Extract reads containing sample barcodes using a barcode probe
Input: 
Raw sequencing reads (FASTQ format)
barcode.probe.fasta

Output: 
ExtractedBC.fasta

Command: 
python UMIC-seq.py UMIextract Plant1_Leaf_UMIC_seq.fq \
  --probe barcode.probe-Plant1.fasta \
  --umi_loc up \
  --umi_len 10 \
  --output Plant1_Leaf_UMIC_ExtractedBC_10bp.fasta &

### Step1.2: Demultiplex reads by sample barcode

In this step, reads in Plant1_Leaf_UMIC_ExtractedBC_10bp.fasta are split according to sample barcodes.
Reads with different barcodes are separated and written to individual FASTA files.

Input: 
barcodes-Plant1.fasta
Plant1_Leaf_UMIC_ExtractedBC_10bp.fasta

output: 
A new directory containing demultiplexed FASTA files, with one FASTA file per sample barcode.

Command: 
nohup python UMIC-seq_helper.py demultiplex \
--barcodes barcodes-Plant1.fasta \
--input Plant1_Leaf_UMIC_ExtractedBC_10bp.fasta \
--output ../Demultiplex_fa/Demultiplex &

### Step1.3: Convert demultiplexed FASTA files to sample-specific FASTQ files

In this step, reads are extracted from Demultiplex_BCXX.fasta and converted into sample-specific FASTQ files based on sample barcodes.

Input: 
Plant1_Leaf_UMIC_seq.fq

output: 
A new directory containing demultiplexed FASTQ files, with one FASTQ file per sample barcode.

Command: 
python UMI-seq-Step1.3.py
See the script UMI-seq-Step1.3.py for implementation details.

## 2. Extraction of UMIs for Each Sample
In this step, Unique Molecular Identifiers (UMIs) are extracted from the sequencing data of each sample.
UMIs are used to distinguish PCR duplicates from true biological variants in downstream analyses.

Input: 
Demultiplexed FASTQ files (../BC_split_fq/BCXX.fq) generated in Step 1.3

Output: 
A new directory containing FASTA files with extracted UMIs, with one FASTA file per sample.

Command:
python UMIC-seq.py UMIextract \
  --input ../BC_split_fq/BCXX.fq \
  --probe UMI_probe-Plant1.fasta \
  --umi_loc down \
  --umi_len 40 \
  --output ../Extract_UMI/

## 3.Full_Clustering

In this step, full UMI clustering is performed.
Reads with identical or highly similar UMIs are grouped together based on the thresholds defined.

Input:  
Extracted UMI FASTA files (Extract_UMI/ExtractedUMIs_BCXX.fasta) generated in Step 2

Output: 
A new directory containing clustered UMIs and associated reads, organized by sample

Command:
python UMIC-seq.py clusterfull \
  --input Extract_UMI/ExtractedUMIs_BCXX.fasta \
  --reads BC_split_fq/BCXX.fq \
  --aln_thresh 50 \
  --size_thresh 10 \
  --output ../UMIclusterfull/ \
  --stop_thresh 0


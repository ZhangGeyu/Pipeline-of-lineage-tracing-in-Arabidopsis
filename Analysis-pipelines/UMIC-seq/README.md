This directory contains the processing pipeline used for UMIC-seq data. This pipeline captures bona fide somatic mutations tracked by UMIs across Arabidopsis samples.
# Pipeline Steps
## 1. **Extraction of Sample Barcodes (Sample_BC)
  This step involves extracting sample barcodes from the raw sequencing data, which are used for identifying and grouping reads according to different samples.
  script: extract_sample_bc.py
  Step1: ### Extraction of reads with sample barcodes by the barcode probe
  Input: Raw sequencing reads (FASTQ format)  and Barcode_Probe.fasta    （放data里）
  Output: ExtractedBC.fasta

Step2: ### Split the Sequence_Data_file by sample barcode 
  Input:  sample_bc.txt  and  ExtractedBC.fasta in step1
  output:  Output_Reads_With_Different_Sample_BC

Step3:### Extract reads from Sequence_Data_file to generate .fastq file for each sample 
  Input:
  output: 
  Refer extract_sample_bc.py script for more details.

## 2. Extraction of UMIs for Each Sample
   In this step, Unique Molecular Identifiers (UMIs) are extracted from each sample’s sequencing data. UMIs help distinguish between PCR duplicates and true biological variants.
   Script: extract_umi.py
   Input: Demultiplexed FASTQ files from Step 1
   Output: data/processed/ExtractedUMIs.fasta (FASTA file containing extracted UMIs for each sample)
   Refer extract_umi.py script for more details.

## 3.Full_Clustering
This step performs the full UMI clustering, where reads with identical or similar UMIs are grouped together based on the thresholds defined.
Script: full_clustering.py
input:  ExtractedUMIs.fasta
output:  
parameters: --aln_thresh 50 --size_thresh 10 --stop_thresh 0
Refer Full_Clustering.py script for more details.

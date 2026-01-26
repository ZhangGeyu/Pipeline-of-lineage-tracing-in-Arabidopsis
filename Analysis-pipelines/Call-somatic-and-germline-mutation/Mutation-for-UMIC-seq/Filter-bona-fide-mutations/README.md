# UMIC-seq Mutation Calling Pipeline
This directory contains the mutation calling pipeline for UMIC-seq experiments.

The pipeline extracts mutations from each read and identifies bona fide somatic mutations for each readout (UMI).

The workflow below uses the Plant1 Leaf UMIC-seq library as an example.

The corresponding input files can be found in: Pipeline-of-lineage-tracing-in-Arabidopsis/Data/UMIC-seq

## Step 1. Merge cluster_XX.fasta Files Produced by Clusterfull
This step merges the clustered FASTA files to facilitate the subsequent mapping step.
#### Input: 
File path of the Clusterfull results (../UMIclusterfull/)
#### Output: 
BC_UMI_merge.fasta
#### Command
python UMIC-seq-CallMut-Step1.py

See the script UMIC-seq-CallMut-Step1.py for implementation details.

## Step 2. Map Reads to Reference and Sort BAM
#### Input: 
BC_UMI_merge.fasta
#### Output: 
BC_UMI_merge.fasta.sorted.bam
#### Command:
minimap2 -ax map-ont -t 50 \
reference-Plant1.fa \
../UMIclusterfull/BC_UMI_merge.fasta \
-o BC_UMI_merge.fasta.bam &

samtools sort BC_UMI_merge.fasta.bam \
-o BC_UMI_merge.fasta.sorted.bam &

## Step 3. Split BAM by Readout
This step splits the sorted BAM so that each readout has its own BAM file.
#### Input: 
BC_UMI_merge.fasta.sorted.bam
#### Output: 
A new directory containing BAM files, one BAM per readout
#### Command:
python UMIC-seq-CallMut-Step3.py

See the script UMIC-seq-CallMut-Step3.py for implementation details.

## Step 4. Call Mutations Using samtools mpileup
#### Input: 
../BAM_and_mpileup_split/BCXX_XX.bam
#### Output: 
../BAM_and_mpileup_split/BCXX_XX.mpileup
#### Command:
samtools mpileup --max-depth 0 --output-BP \
--reference reference-Plant1.fa \
../BAM_and_mpileup_split/BCXX_XX.bam \
-o ../BAM_and_mpileup_split/BCXX_XX.mpileup

## Step 5. Extract Mutations from mpileup Files
This step reformats mutation information from the mpileup file into a dataframe.
#### Input:
SampleBC_SampleName-Plant1.txt
../BAM_and_mpileup_split/
#### Output:
ConsensusSequence_CallSNP_Raw.txt (all mutations from all samples)
#### Command:
python UMIC-seq-CallMut-Step5.py

See the script UMIC-seq-CallMut-Step5.py for implementation details.

## Step 6. Filter Bona Fide Mutations
This step filters sequencing errors and extracts bona fide mutations.
#### Input:
ConsensusSequence_CallSNP_Raw.txt
#### Output:
ConsensusSequence_CallSNP.txt
#### Command:
python UMIC-seq-CallMut-Step5.py
See the script UMIC-seq-CallMut-Step5.py for implementation details.

##### Important: Since one plant may have multiple UMIC-seq libraries, merge all libraries for each plant.
##### The merged file should be named like: Parental_CallSNP_Plant1.txt

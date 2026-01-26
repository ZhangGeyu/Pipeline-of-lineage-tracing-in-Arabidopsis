This directory contains the mutation calling pipeline for UMIC-seq experiments.
The pipeline is designed to extract mutations in each reads, and identify bona fide somatic mutations for each readouts(UMIs).

The workflow below uses the Plant1 Leaf UMIC-seq library as an example.
The corresponding input files can be found in: Pipeline-of-lineage-tracing-in-Arabidopsis/Data/UMIC-seq

## Step 1. Merge the cluster_XX.fasta file produced by Clusterfull
This step is used to make the following mapping step easier
#### Input: 
File path of the Clusterfull results (../UMIclusterfull/)
#### Output: 
BC_UMI_merge.fasta
#### Command
python UMIC-seq-CallMut-Step1.py
See the script UMIC-seq-CallMut-Step1.py for implementation details.

## Step 2. Mapping the readouts to reference and sort bam
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

## Step 3. split the BC_UMI_merge.fasta.sorted.bam and make each readout in a single .bam
#### Input: 
BC_UMI_merge.fasta.sorted.bam
#### Output: 
A new directory containing BAM files, with one BAM file per readout.
#### Command:
python UMIC-seq-CallMut-Step3.py
See the script UMIC-seq-CallMut-Step3.py for implementation details.

## Step 4. Call mutation using samtools mpileup for all BAM
#### Input: 
../BAM_and_mpileup_split/BCXX_XX.bam
#### Output: 
../BAM_and_mpileup_split/BCXX_XX.mpileup
#### Command:
samtools mpileup --max-depth 0 --output-BP \
--reference reference-Plant1.fa \
../BAM_and_mpileup_split/BCXX_XX.bam \
-o ../BAM_and_mpileup_split/BCXX_XX.mpileup

## Step 5. Extract Mutations from mpileup File
This step reformats the mutation information in mpileup file into a dataframe.
#### Input:
SampleBC_SampleName-Plant1.txt
../BAM_and_mpileup_split/
#### Output:
ConsensusSequence_CallSNP_Raw.txt (all of the mutations in all samples)
#### Command:
python UMIC-seq-CallMut-Step5.py
See the script UMIC-seq-CallMut-Step5.py for implementation details.

## Step 6. Filter the bona fide Mutations
This step set a cutoff to Filter the sequencing error, and extract the bona fide Mutations
#### Input:
ConsensusSequence_CallSNP_Raw.txt
#### Output:
ConsensusSequence_CallSNP.txt
#### Command:
python UMIC-seq-CallMut-Step5.py
See the script UMIC-seq-CallMut-Step5.py for implementation details.
#### Note
Because one plant may have several UMIC-seq librarys, please merge all the librarys together for each plant, and the file's name is like "Parental_CallSNP_Plant1.txt"



## 6. Extraction of Readouts for Used Copy
  Extracts specific readouts that correspond to the used readout copy.
  Script: extract_readouts.py
  Input: 
  Output: readouts_used_copy.fasta
  
## 7. Frequency of Mutations in Each Sample
  Calculates the frequency of mutations across each sample.
  Script: compute_mutation_frequency.py
  Input: filtered_mutations.txt
Output: mutation_frequencies.txt

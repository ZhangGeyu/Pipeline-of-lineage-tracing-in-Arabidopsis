# Analysis of Plant 1 Progeny Sequencing (TA-clone)
This section describes the analysis of the Plant 1 progeny sequencing experiment generated using TA cloning followed by Sanger sequencing.

The final output files of this section can be found in: Pipeline-of-lineage-tracing-in-Arabidopsis/Data/Mutation_information

## Step 1. Generate Consensus Sequences from Sanger Sequencing Data
In this step, forward and reverse Sanger sequencing reads are combined to generate a consensus sequence for each TA clone.
#### Input
Raw Sanger sequencing file: Plant1_Progeny_Sanger_seq.fasta
#### Output
Consensus sequences: Plant1_Progeny_merge_seq.fasta
#### Command
python TA-clone-ConsensusSeq.py

## Step 2. Map Consensus Sequences to the Reference and Call Mutations
Consensus sequences are mapped to the reference readout, followed by mutation calling using samtools mpileup
#### Input
Reference sequence: reference-Plant1.fa
Consensus sequences: Plant1_Progeny_merge_seq.fasta
#### Output
mpileup file containing mutation information: Plant1_Progeny_merge_seq.mpileup
#### Command
nohup minimap2 -ax map-ont -t 50 reference-Plant1.fa Plant1_Progeny_merge_seq.fasta -o Plant1_Progeny_merge_seq.fasta.bam &

nohup samtools sort Plant1_Progeny_merge_seq.fasta.bam -o Plant1_Progeny_merge_seq.fasta.sorted.bam &

nohup samtools mpileup --max-depth 0 --output-BP --output-QNAME --reference reference-Plant1.fa Plant1_Progeny_merge_seq.fasta.sorted.bam -o Plant1_Progeny_merge_seq.mpileup &

## Step 3. Extract Mutation Information for Each Clone
This step extracts mutation information for each TA clone from the mpileup file.
#### Input
Plant1_Progeny_merge_seq.mpileup
#### Output
Mutation list for each progeny clone: Progeny_CallSNP_Plant1.txt
#### Command
python TA-clone-Call-mutation.py

## Step 4. Calculate Mutation Frequencies in Each Progeny Sample
This step calculates the mutation frequency across progeny samples.
#### Input
Progeny_CallSNP_Plant1.txt
#### Output
Mutation frequency table: MutFreq_In_ProgenySample_Plant1.txt
#### Command
python TA-clone-mutation-freq.py

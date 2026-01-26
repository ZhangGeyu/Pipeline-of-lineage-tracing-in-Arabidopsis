# Readout Extraction and Mutation Frequency Analysis
This section is used to extract readouts corresponding to the selected copy and to calculate mutation frequencies in each sample for downstream analyses, such as progeny haplotype identification.

The output files generated in this section are frequently used in subsequent analyses.

All results can be found in: Pipeline-of-lineage-tracing-in-Arabidopsis/Data/Mutation_information

## Step 1. Extraction of Readouts for the Selected Copy
In this step, feature mutations are used to extract specific readouts that correspond to the selected copy.

This step is only applied to Plants 2 and 3, for both parental and progeny samples.

Here, the Plant 2 parental sample is used as an example.
#### Input: 
Parental_CallSNP_Plant2.txt
#### Output: 
Parental_CallSNP_Plant2_886copy.txt
#### Command
python UMIC-seq-extract-copy-Plant2.py

See UMIC-seq-extract-copy-Plant2.py for implementation details.

## Step 2. Calculation of Mutation Frequencies in Each Sample
This step calculates the mutation frequency for each mutation across samples.
#### This step is applied to all UMIC-seq libraries, including:
Parental samples of Plants 1, 2, and 3

Progeny samples of Plants 2 and 3

The Plant 2 parental sample is shown here as an example.

#### Input: 
Parental_CallSNP_Plant2_886copy.txt
#### Output: 
MutFreq_In_ParentalSample_Plant2_886copy.txt
#### Command
python UMIC-seq-MutFreq-Plant2.py

See UMIC-seq-MutFreq-Plant2.py for implementation details.

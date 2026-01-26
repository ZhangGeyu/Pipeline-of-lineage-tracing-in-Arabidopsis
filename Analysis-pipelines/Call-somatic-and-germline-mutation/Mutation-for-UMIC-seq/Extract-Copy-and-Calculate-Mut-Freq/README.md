This part is used to extract readouts from the copy we used, and calculate the mutation frequency in each sample for the futher analysis (like progeny haplotype identification)
The output files of this part would be frequently used in the following analysis,the output file could be find at: Pipeline-of-lineage-tracing-in-Arabidopsis/Data/Mutation_information


## Step 1. Extraction of Readouts for Used Copy
Using featured mutation to extract specific readouts that correspond to the used copy.
This step is only used to Plants 2 and 3 for both parental samples and progeny samples, using plant2 parental as the pipeline example.
#### Input: 
Parental_CallSNP_Plant2.txt
#### Output: 
Parental_CallSNP_Plant2_886copy.txt
#### Command
python UMIC-seq-extract-copy-Plant2.py
See UMIC-seq-extract-copy-Plant2.py for implementation details.

## Step 2. Frequency of Mutations in Each Sample
Calculates the frequency of mutations across each sample.
This step is only used to all the UMIC-seq library, including parental samples of Plants 1, 2 and 3, and progeny samples of plant 2 and 3. using plant2 parental as the pipeline example.
#### Input: 
Parental_CallSNP_Plant2_886copy.txt
#### Output: 
MutFreq_In_ParentalSample_Plant2_886copy.txt
#### Command
python UMIC-seq-MutFreq-Plant2.py
See UMIC-seq-MutFreq-Plant2.py for implementation details.

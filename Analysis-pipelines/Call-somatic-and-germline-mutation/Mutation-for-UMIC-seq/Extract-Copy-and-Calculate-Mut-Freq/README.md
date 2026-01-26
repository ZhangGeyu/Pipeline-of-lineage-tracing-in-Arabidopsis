This part is used to extract readouts from the copy we used, and calculate the mutation frequency in each sample for the futher analysis (like progeny haplotype identification)

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

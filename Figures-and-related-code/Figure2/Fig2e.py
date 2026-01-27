# Relationship between the number of unshared mutations for a cauline leaf pair and the number of internodes separating them.
import pandas as pd
import os

# Remove hotspot mutations
Hotspot_Plant1 = ["1092_T_-1C","1114_A_G","1226_C_T","1335_C_T","445_C_T","513_C_G","53_C_-2GA","77_G_A","795_G_T","965_G_A"]

mutation_file = os.path.join('MutFreq_In_ParentalSample_Plant1.txt')
distance_file = os.path.join('sample_distance.txt')
output_file = os.path.join('different_mutation_between_sample_final.txt')

df = pd.read_csv(mutation_file, sep='\t')
df.columns = ['SampleName', 'mut_info', 'mut_info_count'] + list(df.columns[3:])
df['mut_info_count'] = pd.to_numeric(df['mut_info_count'], errors='coerce')
df['SampleName'] = df['SampleName'].astype(str)
df = df[~df['SampleName'].str.contains('RL', case=False, na=False)]

# Use the top 100 mutations for each sample
sample_data = {}
for sample, group in df.groupby('SampleName'):
     # Filter out hotspot mutations
    filtered_group = group[~group['mut_info'].isin(Hotspot_Plant1)]
     # Sort and take top 100 mutations
    top_mutations = group.sort_values('mut_info_count', ascending=False).head(100)
    sample_data[sample] = set(top_mutations['mut_info'].astype(str))

# Read the distance file
distance_map = {}
with open(distance_file, 'r') as f:
    for line in f:
        s1, s2, dist = line.strip().split()
        s1 = s1.replace('_', '-')
        s2 = s2.replace('_', '-')
        key1 = f"{s1}:{s2}"
        key2 = f"{s2}:{s1}"
        distance_map[key1] = dist
        distance_map[key2] = dist


samples = sorted(sample_data.keys())
with open(output_file, 'w', encoding='utf-8') as f:
    for i in range(len(samples)):
        s1 = samples[i]
        for j in range(i + 1, len(samples)):
            s2 = samples[j]
            key = f"{s1}:{s2}"
            if key not in distance_map:
                continue
            diff_count = len(sample_data[s1] ^ sample_data[s2])
            f.write(f"{s1}\t{s2}\t{diff_count}\t{distance_map[key]}\n")

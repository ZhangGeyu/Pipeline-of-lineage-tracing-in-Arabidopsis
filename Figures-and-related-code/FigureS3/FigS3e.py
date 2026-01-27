Hotspot_181copy = ["1007_C_T","1013_C_T","1100_G_A","1186_A_G","1242_G_A","306_T_C","77_C_T","850_C_T","941_C_T"]
HighFreqMut_list_181copy = ['1213_C_+5GTGCT','181_G_A','182_G_A','183_G_A','184_G_A','736_G_A','601_G_A','657_G_A','719_C_T','238_G_A','658_G_A','663_C_T']

import pandas as pd

# Raw file
Offspring_mut_freq = pd.read_table('MutFreq_In_ProgenySample_Plant3_181copy.txt')

Offspring_mut_freq = Offspring_mut_freq[(~Offspring_mut_freq['mut_info'].isin(HighFreqMut_list_181copy))&(~Offspring_mut_freq['mut_info'].isin(Hotspot_181copy))]   # Remove hotspot mutations and copy mutations
Offspring_mut_freq = Offspring_mut_freq[['mut_info','SampleName','mut_freq']]
Offspring_mut_freq = Offspring_mut_freq[Offspring_mut_freq['mut_freq'] >= 0.3]
print(Offspring_mut_freq['SampleName'].nunique())

# Raw file
Offspring_UMI_mut = pd.read_table('Progeny_CallSNP_Plant3_181copy.txt')
Offspring_UMI_mut = Offspring_UMI_mut[(~Offspring_UMI_mut['mut_info'].isin(HighFreqMut_list_181copy))&(~Offspring_UMI_mut['mut_info'].isin(Hotspot_181copy))]   # Remove hotspot mutations and copy mutations
print(Offspring_UMI_mut['SampleName'].nunique())

rows = []

for i in Offspring_mut_freq['SampleName'].unique():
    Offspring_mut_freq_subset = Offspring_mut_freq[Offspring_mut_freq['SampleName'] == i]
    Offspring_UMI_mut_subset = Offspring_UMI_mut[Offspring_UMI_mut['SampleName'] == i]

    mut_list = Offspring_mut_freq_subset['mut_info'].unique()
    Offspring_UMI_mut_subset = Offspring_UMI_mut_subset[Offspring_UMI_mut_subset['mut_info'].isin(mut_list)]
    if Offspring_UMI_mut_subset['mut_info'].nunique() >= 1:
        Offspring_UMI_mut_subset = Offspring_UMI_mut_subset[['SampleName_UMI','mut_info']]
        Offspring_UMI_mut_subset.to_csv((path + 'Sample_Haplotype/SpermAndEgg/' + i + '.txt'), sep = '\t', index = False)

        umi_mutset = (Offspring_UMI_mut_subset.groupby('SampleName_UMI')['mut_info'].apply(lambda x: tuple(sorted(set(x)))).reset_index(name='mut_combo'))
        combo_counts = (umi_mutset.groupby('mut_combo').size().reset_index(name='UMI_count').sort_values('UMI_count', ascending=False)).reset_index(drop = True)
        total_umi = umi_mutset['SampleName_UMI'].nunique()

        top2 = combo_counts.head(2).copy()
        top2['fraction'] = top2['UMI_count'] / total_umi
        result = top2[top2['fraction'] > 0.3]

        for j in range(len(result['mut_combo'])):
            if j == 0:
                mut_combo = result['mut_combo'][j]
                for k in range(len(mut_combo)):
                    rows.append({'mut_info': mut_combo[k],'SampleName': i,'SampleName_UMI': (i + '_1st'),\
                        'haplotype': mut_combo,'UMI_count':result['UMI_count'][j],'fraction':result['fraction'][j]})
            elif j == 1:
                mut_combo = result['mut_combo'][j]
                for k in range(len(mut_combo)):
                    rows.append({'mut_info': mut_combo[k],'SampleName': i,'SampleName_UMI': (i + '_2nd'),\
                        'haplotype': mut_combo,'UMI_count':result['UMI_count'][j],'fraction':result['fraction'][j]})
    #else:
        #print(i)

Offspring_haplotype_df = pd.DataFrame(rows)

# Clean file
Offspring_haplotype_df.to_csv('Offspring_haplotype-Plant3.txt', sep = '\t',index = False)

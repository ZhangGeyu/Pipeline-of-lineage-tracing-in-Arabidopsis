# Plant3

Hotspot_181copy = ["1007_C_T","1013_C_T","1100_G_A","1186_A_G","1242_G_A","306_T_C","77_C_T","850_C_T","941_C_T"]

import pandas as pd

# Raw file
MutFreq_In_ProgenySample = pd.read_table('MutFreq_In_ProgenySample_Plant3_181copy.txt')
MutFreq_In_ProgenySample = MutFreq_In_ProgenySample[~MutFreq_In_ProgenySample['mut_info'].isin(Hotspot_181copy)]

# Raw file
Offspring_haplotype_df = pd.read_table('Offspring_haplotype-Plant3.txt')

to_remove = pd.MultiIndex.from_frame(Offspring_haplotype_df[['SampleName', 'mut_info']])
PostMut_In_ProgenySample = MutFreq_In_ProgenySample[~pd.MultiIndex.from_frame(MutFreq_In_ProgenySample[['SampleName', 'mut_info']]).isin(to_remove)]

PostMut_In_ProgenySample = PostMut_In_ProgenySample[PostMut_In_ProgenySample['mut_info_count'] >=2]

PostMut_In_ProgenySample['ref'] = PostMut_In_ProgenySample['mut_info'].str.split('_').str.get(1).str.upper()
PostMut_In_ProgenySample['alt'] = PostMut_In_ProgenySample['mut_info'].str.split('_').str.get(2).str.upper()

PostMut_In_ProgenySample['Mut_type'] = np.where(PostMut_In_ProgenySample['alt'].str.contains('\+'),'insertion',
    np.where(PostMut_In_ProgenySample['alt'].str.contains('-'),'deletion',
    PostMut_In_ProgenySample['ref'] + '>' + PostMut_In_ProgenySample['alt']))

Mut_Type_Count = PostMut_In_ProgenySample.value_counts(['Mut_type']).reset_index()
Mut_Type_Count.columns = ['Mut_type','count']

# Clean file
with open('MutTypeCount_181copy_Post-germination.txt', 'w') as MutTypeCount:
    MutTypeCount.write('Mut_type\t'+'count\n')
    # Define mutation types
    for i in ['C>T/G>A','C>G/G>C','C>A/G>T','T>A/A>T','T>C/A>G','T>G/A>C']:
        try:
            count_1 = Mut_Type_Count[Mut_Type_Count['Mut_type'] == i.split('/')[0]]['count'].unique()[0]
        except IndexError:
            count_1 = 0   
        try:
            count_2 = Mut_Type_Count[Mut_Type_Count['Mut_type'] == i.split('/')[1]]['count'].unique()[0]
        except IndexError:
            count_2 = 0
        MutTypeCount.write(i + '\t' + str(count_1 + count_2) + '\n')

# Plant2

import numpy as np

Hotspot_886copy = [""]

# Raw file
MutFreq_In_ProgenySample = pd.read_table('MutFreq_In_ProgenySample_Plant2_886copy.txt')

MutFreq_In_ProgenySample = MutFreq_In_ProgenySample[MutFreq_In_ProgenySample['mut_info_count'] >=2]
MutFreq_In_ProgenySample = MutFreq_In_ProgenySample[(MutFreq_In_ProgenySample['mut_freq'] < 0.5)&(~MutFreq_In_ProgenySample['mut_info'].isin(Hotspot_886copy))] # Further filtering: keep records with mutation frequency <0.5

MutFreq_In_ProgenySample['ref'] = MutFreq_In_ProgenySample['mut_info'].str.split('_').str.get(1).str.upper()
MutFreq_In_ProgenySample['alt'] = MutFreq_In_ProgenySample['mut_info'].str.split('_').str.get(2).str.upper()
# Classify mutations by type: insertion, deletion, or base substitution
MutFreq_In_ProgenySample['Mut_type'] = np.where(MutFreq_In_ProgenySample['alt'].str.contains('\+'),'insertion',
    np.where(MutFreq_In_ProgenySample['alt'].str.contains('-'),'deletion',
    MutFreq_In_ProgenySample['ref'] + '>' + MutFreq_In_ProgenySample['alt']))
#  Count occurrences of each mutation type
Mut_Type_Count = MutFreq_In_ProgenySample.value_counts(['Mut_type']).reset_index()
Mut_Type_Count.columns = ['Mut_type','count']

# Clean file
with open('MutTypeCount_886copy_Post-germination.txt', 'w') as MutTypeCount:
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

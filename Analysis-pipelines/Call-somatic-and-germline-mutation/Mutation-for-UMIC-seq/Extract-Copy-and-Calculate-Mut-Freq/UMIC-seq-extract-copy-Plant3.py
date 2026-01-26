import pandas as pd

HighFreqMut_list_181copy = ['1213_C_+5GTGCT','181_G_A','182_G_A','183_G_A','184_G_A','736_G_A','601_G_A','657_G_A','719_C_T','238_G_A','658_G_A','663_C_T']

path = '../' # path of the folder contains Parental_CallSNP_Plant3.txt
copy = '181copy'

Raw_Mutation = pd.read_table((path + 'Parental_CallSNP_Plant3.txt'))

Raw_Mutation_copy = Raw_Mutation[Raw_Mutation['mut_info'].isin(locals()['HighFreqMut_list_' + copy])]
Set_UMI_mutCount = Raw_Mutation_copy.value_counts('SampleName_UMI').reset_index(name = 'count')
overlap_set = Set_UMI_mutCount[Set_UMI_mutCount['count'] == len(locals()['HighFreqMut_list_' + copy])]['SampleName_UMI'].unique()
print(len(overlap_set))

Mut_copy = Raw_Mutation[(Raw_Mutation['SampleName_UMI'].isin(overlap_set))][['pos','ref','alt','SampleName','UMI','SampleName_UMI','mut_info']].reset_index(drop = True)

# Sequences around this deletion cause inaccurate mutation calling
# This site will be excluded from all downstream analyses
Mut_copy = Mut_copy[Mut_copy['mut_info'] != '762_G_-15AGCAACTAATCTAAT']

Mut_copy.to_csv((path + 'Parental_CallSNP_Plant3_' + copy + '.txt'), sep = '\t',index = False)
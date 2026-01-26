import pandas as pd

HighFreqMut_list_886copy = ['1213_C_+5GTGTG','886_G_A','862_G_C','874_G_A','837_G_C','904_G_T','841_C_T','458_T_-20AACAGGGTAATGAGCCGCAC']

path = '../' # path of the folder contains Parental_CallSNP_Plant2.txt
copy = '886copy'

Raw_Mutation = pd.read_table((path + 'Parental_CallSNP_Plant2.txt'))

Raw_Mutation_copy = Raw_Mutation[Raw_Mutation['mut_info'].isin(locals()['HighFreqMut_list_' + copy])]
Set_UMI_mutCount = Raw_Mutation_copy.value_counts('SampleName_UMI').reset_index(name = 'count')
overlap_set = Set_UMI_mutCount[Set_UMI_mutCount['count'] == len(locals()['HighFreqMut_list_' + copy])]['SampleName_UMI'].unique()
print(len(overlap_set))

Mut_copy = Raw_Mutation[(Raw_Mutation['SampleName_UMI'].isin(overlap_set))][['pos','ref','alt','SampleName','UMI','SampleName_UMI','mut_info']].reset_index(drop = True)
Mut_copy.to_csv((path + 'Parental_CallSNP_Plant2' + copy + '.txt'), sep = '\t',index = False)
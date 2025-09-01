## Especially for Plant2 and Plant3

# Plant 2
HighFreqMut_list = ['1213_C_+5GTGTG','886_G_A','862_G_C','874_G_A','837_G_C','904_G_T','841_C_T','458_T_-20AACAGGGTAATGAGCCGCAC']
# Plant 3
HighFreqMut_list = ['1213_C_+5GTGCT','181_G_A','182_G_A','183_G_A','184_G_A','736_G_A','601_G_A','657_G_A','719_C_T','238_G_A','658_G_A','663_C_T']

Clean_Mutation_file = ''
Mutation_file_For_Copy_Used = ''

Raw_Mutation = pd.read_table(Clean_Mutation_file)
Raw_Mutation['mut_info'] = Raw_Mutation['pos'].astype(str) + '_' + Raw_Mutation['ref'] + '_' + Raw_Mutation['alt']

Raw_Mutation_copy = Raw_Mutation[Raw_Mutation['mut_info'].isin(HighFreqMut_list)]
Set_UMI_mutCount_copy = Raw_Mutation_copy.value_counts('SampleName_UMI').reset_index()
overlap_set_copy = Set_UMI_mutCount_copy[Set_UMI_mutCount_copy[0] == len(HighFreqMut_list)]['SampleName_UMI'].unique()
print(copy + ' UMI count: ' + str(len(overlap_set_copy)))

Mut_copy = Raw_Mutation[Raw_Mutation['SampleName_UMI'].isin(overlap_set_copy)]
Mut_copy.to_csv(Mutation_file_For_Copy_Used, sep = '\t',index = False)
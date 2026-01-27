import pandas as pd

# Plant3

Hotspot_Plant3 = ["1007_C_T","1013_C_T","1100_G_A","1186_A_G","1242_G_A","306_T_C","77_C_T","850_C_T","941_C_T"]
HighFreqMut_list_181copy = ['1213_C_+5GTGCT','181_G_A','182_G_A','183_G_A','184_G_A','736_G_A','601_G_A','657_G_A','719_C_T','238_G_A','658_G_A','663_C_T']

# Raw file
Germline_mutation = pd.read_table('Offspring_haplotype.txt')
Germline_mutation = Germline_mutation[(~Germline_mutation['mut_info'].isin(Hotspot_Plant3))&(~Germline_mutation['mut_info'].isin(HighFreqMut_list_181copy))][['mut_info']]
Germline_mutation['group'] = 'Progeny'
Germline_mut_list = Germline_mutation['mut_info'].unique()

# Raw file
Parental_mutation = pd.read_table('Parental_CallSNP_Plant3_181copy.txt')
Parental_mutation = Parental_mutation[(~Parental_mutation['mut_info'].isin(Hotspot_Plant3))&(~Parental_mutation['mut_info'].isin(HighFreqMut_list_181copy))]

Parental_mutation = Parental_mutation[~Parental_mutation['SampleName'].str.contains('RL')] # remove rosette leaves
Parental_mutation['Branch'] = Parental_mutation['SampleName'].str.get(1)
Parental_Branch_list = Parental_mutation['Branch'].unique()
print(Parental_Branch_list)

AllBranch_combinations = [sorted(list(combo)) for r in range(1, len(Parental_Branch_list)+1) for combo in combinations(Parental_Branch_list, r)]
AllBranch_combinations = ['_'.join(str(x) for x in sublist) for sublist in AllBranch_combinations]

Parental_mut_list = Parental_mutation['mut_info'].unique()

PassToOffspring_df = pd.DataFrame({'Parent': AllBranch_combinations,\
    'Parent_Mut': [0] * len(AllBranch_combinations),\
        'PassToOffspring': [0] * len(AllBranch_combinations)})
PassToOffspring_df = PassToOffspring_df.set_index('Parent')

for mut in Parental_mut_list:
    Parental_mutation_subset = Parental_mutation[Parental_mutation['mut_info'] == mut]
    subset_Branch = sorted(list(Parental_mutation_subset['Branch'].unique()))
    subset_Branch = '_'.join(str(x) for x in subset_Branch)
    #print(subset_Branch)
    PassToOffspring_df['Parent_Mut'][subset_Branch] = PassToOffspring_df['Parent_Mut'][subset_Branch] + 1
    if mut in Germline_mut_list:
        PassToOffspring_df['PassToOffspring'][subset_Branch] = PassToOffspring_df['PassToOffspring'][subset_Branch] + 1
    else:
        continue

PassToOffspring_df = PassToOffspring_df.reset_index()
PassToOffspring_df['Fraction'] = PassToOffspring_df['PassToOffspring'] / PassToOffspring_df['Parent_Mut']

PassToOffspring_df['BranchCount'] = PassToOffspring_df['Parent'].str.split('_').apply(len).astype(str) + '_Branch_Shared'

Parent_df = PassToOffspring_df.groupby('BranchCount')['Parent_Mut'].sum().reset_index()
Offspring_df = PassToOffspring_df.groupby('BranchCount')['PassToOffspring'].sum().reset_index()

PassToOffspring_df = pd.merge(Parent_df,Offspring_df)
PassToOffspring_df['Fraction'] = PassToOffspring_df['PassToOffspring'] / PassToOffspring_df['Parent_Mut']

# Clean file
PassToOffspring_df.to_csv('Mut_PassToOffspring_BranchShared-Plant3.txt', sep = '\t', index = False)
PassToOffspring_df

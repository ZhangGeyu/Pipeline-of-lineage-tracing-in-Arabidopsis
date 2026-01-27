# Plant 1

from itertools import combinations

### Parental Pass to Offspring

Hotspot_Plant1 = ["1092_T_-1C","1114_A_G","1226_C_T","1335_C_T","445_C_T","513_C_G","53_C_-2GA","77_G_A","795_G_T","965_G_A"]

# Raw file
Germline_mutation = pd.read_table('MutFreq_In_ProgenySample_Plant1.txt')
Germline_mutation = Germline_mutation[(~Germline_mutation['mut_info'].isin(Hotspot_Plant1))&(Germline_mutation['mut_freq'] >= 0.5)].drop_duplicates('mut_info')[['mut_info']]
Germline_mutation['group'] = 'Progeny'
Germline_mut_list = Germline_mutation['mut_info'].unique()

# Raw file
Parental_mutation = pd.read_table('Parental_CallSNP_Plant1.txt')
Parental_mutation = Parental_mutation[~Parental_mutation['mut_info'].isin(Hotspot_Plant1)]

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
PassToOffspring_df.to_csv('Mut_PassToOffspring_BranchShared-Plant1.txt', sep = '\t', index = False)
PassToOffspring_df

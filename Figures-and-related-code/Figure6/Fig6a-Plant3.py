import pandas as pd

# Plant3

Hotspot_Plant3 = ["1007_C_T","1013_C_T","1100_G_A","1186_A_G","1242_G_A","306_T_C","77_C_T","850_C_T","941_C_T"]
HighFreqMut_list_181copy = ['1213_C_+5GTGCT','181_G_A','182_G_A','183_G_A','184_G_A','736_G_A','601_G_A','657_G_A','719_C_T','238_G_A','658_G_A','663_C_T']

# Raw file
Germline_mutation = pd.read_table('Offspring_haplotype-Plant3.txt')
Germline_mutation = Germline_mutation[(~Germline_mutation['mut_info'].isin(Hotspot_Plant3))&(~Germline_mutation['mut_info'].isin(HighFreqMut_list_181copy))][['mut_info']]
Germline_mutation['group'] = 'Progeny'

# Raw file
Parental_mutation = pd.read_table('Parental_CallSNP_Plant3_181copy.txt')
Parental_mutation = Parental_mutation[(~Parental_mutation['mut_info'].isin(Hotspot_Plant3))&(~Parental_mutation['mut_info'].isin(HighFreqMut_list_181copy))]

Parental_mutation_L = Parental_mutation[Parental_mutation['SampleName'].str.contains('CL')].drop_duplicates('mut_info')[['mut_info']]
Parental_mutation_L['group'] = 'Cauline_Leaves'
Parental_mutation_R = Parental_mutation[Parental_mutation['SampleName'].str.contains('RL')].drop_duplicates('mut_info')[['mut_info']]
Parental_mutation_R['group'] = 'Rosette_Leaves'

All_Mutation = pd.concat([Germline_mutation,Parental_mutation_L,Parental_mutation_R])

# Clean file
All_Mutation.to_csv((path + 'Germline/Consensus_seq/MutationOverview/Parental_Progeny_MutOverlap-Plant3.txt'), sep = '\t', index = False)

import pandas as pd

# Plant2

Hotspot_Plant2 = [""]
HighFreqMut_list_886copy = ['1213_C_+5GTGTG','886_G_A','862_G_C','874_G_A','837_G_C','904_G_T','841_C_T','458_T_-20AACAGGGTAATGAGCCGCAC']

# Raw file
Germline_mutation = pd.read_table('MutFreq_In_ProgenySample_Plant2_886copy.txt')
Germline_mutation = Germline_mutation[(Germline_mutation['mut_freq'] >= 0.5)&\
    (~Germline_mutation['mut_info'].isin(Hotspot_Plant2))&(~Germline_mutation['mut_info'].isin(HighFreqMut_list_886copy))][['mut_info']]
Germline_mutation['group'] = 'Progeny'

# Raw file
Parental_mutation = pd.read_table('Parental_CallSNP_Plant2_886copy.txt')
Parental_mutation = Parental_mutation[(~Parental_mutation['mut_info'].isin(Hotspot_Plant2))&(~Parental_mutation['mut_info'].isin(HighFreqMut_list_886copy))]

Parental_mutation_L = Parental_mutation[Parental_mutation['SampleName'].str.contains('CL')].drop_duplicates('mut_info')[['mut_info']]
Parental_mutation_L['group'] = 'Cauline_Leaves'
Parental_mutation_R = Parental_mutation[Parental_mutation['SampleName'].str.contains('RL')].drop_duplicates('mut_info')[['mut_info']]
Parental_mutation_R['group'] = 'Rosette_Leaves'

All_Mutation = pd.concat([Germline_mutation,Parental_mutation_L,Parental_mutation_R])

# Clean file
All_Mutation.to_csv('Parental_Progeny_MutOverlap-Plant2.txt', sep = '\t', index = False)

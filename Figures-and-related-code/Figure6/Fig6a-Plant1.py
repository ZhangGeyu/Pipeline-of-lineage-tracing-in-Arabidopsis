import pandas as pd

# Plant1

Hotspot_Plant1 = ["1092_T_-1C","1114_A_G","1226_C_T","1335_C_T","445_C_T","513_C_G","53_C_-2GA","77_G_A","795_G_T","965_G_A"]

# Raw file
Germline_mutation = pd.read_table('MutFreq_In_ProgenySample_Plant1.txt')
Germline_mutation = Germline_mutation[(~Germline_mutation['mut_info'].isin(Hotspot_Plant1))&(Germline_mutation['mut_freq'] >= 0.5)].drop_duplicates('mut_info')[['mut_info']]
Germline_mutation['group'] = 'Progeny'

# Raw file
Parental_mutation = pd.read_table('Parental_CallSNP_Plant1.txt')
Parental_mutation = Parental_mutation[~Parental_mutation['mut_info'].isin(Hotspot_Plant1)]

Parental_mutation_L = Parental_mutation[Parental_mutation['SampleName'].str.contains('CL')].drop_duplicates('mut_info')[['mut_info']]
Parental_mutation_L['group'] = 'Cauline_Leaves'
Parental_mutation_R = Parental_mutation[Parental_mutation['SampleName'].str.contains('RL')].drop_duplicates('mut_info')[['mut_info']]
Parental_mutation_R['group'] = 'Rosette_Leaves'

All_Mutation = pd.concat([Germline_mutation,Parental_mutation_L,Parental_mutation_R])

# Clean file
All_Mutation.to_csv('Parental_Progeny_MutOverlap-Plant1.txt', sep = '\t', index = False)

import pandas as pd

Hotspot_Plant1 = ["1092_T_-1C","1114_A_G","1226_C_T","1335_C_T","445_C_T","513_C_G","53_C_-2GA","77_G_A","795_G_T","965_G_A"]

# Raw file for mutation information
Progeny_mutation = pd.read_table('Progeny_CallSNP_Plant1.txt')
Progeny_mutation = Progeny_mutation[~Progeny_mutation['mut_info'].isin(Hotspot_Plant1)] # remove hotspot mutations

# Raw file contains position of the target site
TargetSite = pd.read_table(path + 'TargetSite.txt')
TargetSite_List = list(TargetSite['position'].unique())
# find nearest target site
def find_nearest_value(row):
    return min(TargetSite_List, key=lambda x: abs(row['pos'] - x))

Progeny_mutation['TargetSite'] = Progeny_mutation.apply(find_nearest_value, axis=1)
Progeny_mutation['Distance'] = Progeny_mutation['pos'] - Progeny_mutation['TargetSite']

TargetSite.columns = ['Target','TargetSite','direction']
Progeny_mutation = pd.merge(Progeny_mutation,TargetSite,on = 'TargetSite')

for i in range(len(Progeny_mutation['pos'])):
    if Progeny_mutation['direction'][i] == '<=':
        Progeny_mutation['Distance'][i] = -Progeny_mutation['Distance'][i]

# Clean file
Progeny_mutation[['SampleName_UMI','mut_info','Distance']].to_csv('DistanceFromTargetSite.txt', sep = '\t', index = False)

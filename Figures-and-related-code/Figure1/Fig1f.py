import pandas as pd

Hotspot_Plant1 = ["1092_T_-1C","1114_A_G","1226_C_T","1335_C_T","445_C_T","513_C_G","53_C_-2GA","77_G_A","795_G_T","965_G_A"]

# Raw file for Parental_mutation
Parental_mutation = pd.read_table('Parental_CallSNP_Plant1.txt')
Parental_mutation = Parental_mutation[~Parental_mutation['mut_info'].isin(Hotspot_Plant1)] # remove hotspot mutations

# Raw file contains the target position
TargetSite = pd.read_table('TargetSite.txt')
TargetSite_List = list(TargetSite['position'].unique())
# find nearest target site
def find_nearest_value(row):
    return min(TargetSite_List, key=lambda x: abs(row['pos'] - x))

Parental_mutation['TargetSite'] = Parental_mutation.apply(find_nearest_value, axis=1)
Parental_mutation['Distance'] = Parental_mutation['pos'] - Parental_mutation['TargetSite']

TargetSite.columns = ['Target','TargetSite','direction']
Parental_mutation = pd.merge(Parental_mutation,TargetSite,on = 'TargetSite')

for i in range(len(Parental_mutation['pos'])):
    if Parental_mutation['direction'][i] == '<=':
        Parental_mutation['Distance'][i] = -Parental_mutation['Distance'][i]

# Clean file
Parental_mutation[['SampleName_UMI','mut_info','Distance']].to_csv('DistanceFromTargetSite.txt', sep = '\t', index = False)

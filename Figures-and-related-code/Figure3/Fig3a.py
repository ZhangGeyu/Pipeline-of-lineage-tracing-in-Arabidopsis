# Observed number of germline mutations shared among branches of Plant 1.
# Remove hotspot mutation
Hotspot_Plant1 = ["1092_T_-1C","1114_A_G","1226_C_T","1335_C_T","445_C_T","513_C_G","53_C_-2GA","77_G_A","795_G_T","965_G_A"]

# Raw file
Germline_mutation = pd.read_table('MutFreq_In_ProgenySample_Plant1.txt')
Germline_mutation = Germline_mutation[(~Germline_mutation['mut_info'].isin(Hotspot_Plant1))&(Germline_mutation['mut_freq'] >= 0.5)]

# Annotate mutations by branch.
Germline_Branch_Mut = Germline_mutation.copy()
Germline_Branch_Mut['Branch'] = 'Branch_' + Germline_Branch_Mut['SampleName'].str.get(1)
Germline_Branch_Mut['Branch_mut'] =  Germline_Branch_Mut['Branch'] + '-' + Germline_Branch_Mut['mut_info']
Germline_Branch_Mut = Germline_Branch_Mut.drop_duplicates('Branch_mut')

# Clean file
Germline_Branch_Mut[['mut_info','Branch']].to_csv('Progeny-BranchMutVenn.txt', sep = '\t', index = False)

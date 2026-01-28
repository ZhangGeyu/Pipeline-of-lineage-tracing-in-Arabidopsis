# Plant3

import pandas as pd

Hotspot_Plant3 = ["1007_C_T","1013_C_T","1100_G_A","1186_A_G","1242_G_A","306_T_C","77_C_T","850_C_T","941_C_T"]
HighFreqMut_list_181copy = ['1213_C_+5GTGCT','181_G_A','182_G_A','183_G_A','184_G_A','736_G_A','601_G_A','657_G_A','719_C_T','238_G_A','658_G_A','663_C_T']
# Raw file
Germline_mutation = pd.read_table('Offspring_haplotype-Plant3.txt')
Germline_mutation = Germline_mutation[(~Germline_mutation['mut_info'].isin(Hotspot_Plant3))&(~Germline_mutation['mut_info'].isin(HighFreqMut_list_181copy))][['mut_info','SampleName','SampleName_UMI']]
# Extract position, reference base, and alternative base from mutation info
Germline_mutation['pos'] = Germline_mutation['mut_info'].str.split('_').str.get(0).astype(int)
Germline_mutation['ref'] = Germline_mutation['mut_info'].str.split('_').str.get(1)
Germline_mutation['alt'] = Germline_mutation['mut_info'].str.split('_').str.get(2)
Germline_mutation['Ref_Alt'] = Germline_mutation['ref'] + '_' + Germline_mutation['alt']

# Raw file: parental SNP data
Parental_mutation = pd.read_table('Parental_CallSNP_Plant3_181copy.txt')
Parental_mutation = Parental_mutation[(~Parental_mutation['mut_info'].isin(Hotspot_Plant3))&\
    (~Parental_mutation['mut_info'].isin(HighFreqMut_list_181copy))][['pos','ref','alt','mut_info','SampleName','SampleName_UMI']]

Parental_mutation['Ref_Alt'] = Parental_mutation['ref'] + '_' + Parental_mutation['alt']

# State for each mutation

Merge_mut_df = pd.concat([Parental_mutation,Germline_mutation]).reset_index(drop = True)  # Merge parental and germline mutation data
Merge_mut_df['ref'] = Merge_mut_df['ref'].str.upper()
Merge_mut_df['alt'] = Merge_mut_df['alt'].str.upper()
Merge_mut_df = Merge_mut_df[Merge_mut_df['alt'].isin(['A','C','G','T'])] # Only SNP

# Give each ref_alt a state number
Ref_Alt_list = list(Merge_mut_df['Ref_Alt'].unique())
Ref_Alt_dict = dict(zip(Ref_Alt_list, range(1, len(Ref_Alt_list) + 1)))
Ref_Alt_df = pd.DataFrame.from_dict(Ref_Alt_dict, orient='index').reset_index()
Ref_Alt_df.columns = ['Ref_Alt','state']

Parental_mutation = pd.merge(Parental_mutation, Ref_Alt_df, on = 'Ref_Alt')
Germline_mutation = pd.merge(Germline_mutation, Ref_Alt_df, on = 'Ref_Alt')
Merge_mut_df = pd.merge(Merge_mut_df, Ref_Alt_df, on = 'Ref_Alt')

### generate mutation_prior.csv

Total_UMI_count = Merge_mut_df['SampleName_UMI'].nunique()
pos_dict = {i: f'c{i-1}' for i in range(1, 1354)}
pos_df = pd.DataFrame.from_dict(pos_dict, orient='index').reset_index()
pos_df.columns = ['pos','character']

mutation_prior_df = pd.DataFrame(columns=['character','state','probability'])

for pos in Merge_mut_df['pos'].unique():
    Merge_mut_df_pos = Merge_mut_df[Merge_mut_df['pos'] == pos]
    C_character = pos_dict[pos]

    state_count = Merge_mut_df_pos.value_counts('state').reset_index()
    state_count.columns = ['state','count']
    state_count['probability'] = state_count['count'] / Total_UMI_count
    state_count['character'] = C_character
    state_count = state_count[['character','state','probability']]
    mutation_prior_df = pd.concat([mutation_prior_df, state_count]).reset_index(drop = True)

## using the tree constructed by FastTree
from skbio import TreeNode
tree = TreeNode.read("ConsensusSeq_Tree_Plant3.nwk")
Parental_UMI_list = [n.name.rsplit(' ', 1)[0].replace(' ', '-') + '_' + n.name.rsplit(' ', 1)[1] if ' ' in n.name else n.name for n in tree.tips()]
Parental_UMI_list = [s for s in Parental_UMI_list if 'P' not in s]

Parental_mutation = Parental_mutation[Parental_mutation['SampleName_UMI'].isin(Parental_UMI_list)].reset_index(drop = True)
print('Parental UMI count: ' + str(Parental_mutation['SampleName_UMI'].nunique()))

Merge_mut_df = pd.concat([Parental_mutation, Germline_mutation]).reset_index(drop = True)
Merge_mut_df = pd.merge(Merge_mut_df, pos_df, on = 'pos')

P_O_mut_Matrix = Merge_mut_df.pivot_table(index='SampleName_UMI', columns='character', values='state', fill_value=0)
P_O_mut_Matrix.index.name = None
P_O_mut_Matrix.columns.name = None
P_O_mut_Matrix.loc['reference'] = 0
Merge_mut_df['CS'] = Merge_mut_df['character'] + '_' + Merge_mut_df['state'].astype(str)

mutation_prior_df = mutation_prior_df.merge(Merge_mut_df.drop_duplicates('CS')[['state','character']], on=['state', 'character'], how='inner').sort_values(['character','state'])
character_50_list = list(mutation_prior_df['character'].unique())[0:600]

P_O_mut_Matrix = P_O_mut_Matrix.loc[:, P_O_mut_Matrix.columns.isin(character_50_list)]
P_O_mut_Matrix = P_O_mut_Matrix.astype(int)
P_O_mut_Matrix.index = P_O_mut_Matrix.index.str.replace("_", "-", regex=False)
P_O_mut_Matrix.to_csv('character_matrix.csv')

mutation_prior_df = mutation_prior_df[mutation_prior_df['character'].isin(character_50_list)]
mutation_prior_df.to_csv('mutation_prior.csv', index = False)

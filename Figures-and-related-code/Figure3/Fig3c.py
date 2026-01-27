import pandas as pd

# Example of mutations shared among progenies from different plant branches.
# Remove hotspot mutation

# Linkage mutations list
MutSet_list = ['379_C_T','439_C_T','630_C_T','981_C_T']
str_linkage = str(MutSet_list).replace('_','').replace('[','').replace(']','').replace("'",'').replace(',','_').replace(' ','')

# Raw file of parental mutations
Parental_Mut = pd.read_table('Parental_CallSNP_Plant1.txt')

Sample_UMI_count = Parental_Mut.drop_duplicates('SampleName_UMI').value_counts(['SampleName']).reset_index()
Sample_UMI_count.columns = ['SampleName','UMI_count']

# Raw file of germline mutations
Offspring_Mut = pd.read_table('MutFreq_In_ProgenySample_Plant1.txt')
Offspring_Mut = Offspring_Mut[Offspring_Mut['mut_freq'] >= 0.5]

# Clean file  
# Co-occurrence frequency of linked mutations across different samples.
with open(('Example-' + str_linkage + '.txt'),'w') as MutSet_freq:
    MutSet_freq.write('sample\t'+'mut_info\t'+'Linkage_UMI_count\t'+'UMI_count\t'+'mut_freq\n')
    
    linkage_df = Parental_Mut[Parental_Mut['mut_info'].isin(MutSet_list)]
    UMI_MutCount = linkage_df.value_counts('SampleName_UMI').reset_index(name = 'count')
    
    linkage_overlap_set = list(UMI_MutCount[UMI_MutCount['count'] == len(MutSet_list)]['SampleName_UMI'].unique())

    linkage_UMI_df = pd.DataFrame(linkage_overlap_set, columns=['SampleName_UMI'])
    linkage_UMI_df['SampleName'] = linkage_UMI_df['SampleName_UMI'].str.split('_').str.get(0)
    linkage_SampleName_count = linkage_UMI_df.value_counts('SampleName').reset_index()
    linkage_SampleName_count.columns = ['SampleName','Linkage_UMI_count']

    Offspring_linkage_df = Offspring_Mut[Offspring_Mut['mut_info'].isin(MutSet_list)]
    Offspring_linkage_count = Offspring_linkage_df.value_counts('SampleName').reset_index(name = 'count')
    Offspring_SampleName_list = Offspring_linkage_count[Offspring_linkage_count['count'] == len(MutSet_list)]['SampleName'].unique()
    Offspring_linkage_UMI_df = pd.DataFrame(Offspring_SampleName_list, columns=['SampleName'])
    Offspring_linkage_UMI_df['Linkage_UMI_count'] = 1

    linkage_SampleName_count = pd.merge(linkage_SampleName_count, Sample_UMI_count, how = 'outer')

    Offspring_linkage_UMI_df['UMI_count'] = 1

    linkage_SampleName_count = pd.concat([linkage_SampleName_count, Offspring_linkage_UMI_df]).reset_index(drop = True)

    linkage_SampleName_count = linkage_SampleName_count.fillna(0)
    linkage_SampleName_count['UMI_fraction'] = linkage_SampleName_count['Linkage_UMI_count'] / linkage_SampleName_count['UMI_count']
    linkage_SampleName_count = linkage_SampleName_count.sort_values('UMI_fraction', ascending = False).reset_index(drop = True)
        
    linkage_Parental = linkage_SampleName_count[linkage_SampleName_count['UMI_fraction'] != 1]
    for i in range(len(linkage_SampleName_count['SampleName'])):
        MutSet_freq.write(linkage_SampleName_count['SampleName'][i] + '\t' + str_linkage + '\t' + str(linkage_SampleName_count['Linkage_UMI_count'][i]) + '\t' + \
            str(linkage_SampleName_count['UMI_count'][i]) + '\t' + str(linkage_SampleName_count['UMI_fraction'][i]) + '\n')

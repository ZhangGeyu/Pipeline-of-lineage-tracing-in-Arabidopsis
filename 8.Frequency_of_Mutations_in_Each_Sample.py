import pandas as pd

Mutation_file_For_Copy_Used = ''
Mutation_Frequency_file = ''

Mut_total = pd.read_table(Mutation_file_For_Copy_Used)

Mut_total['ref'] = Mut_total['ref'].str.upper()
Mut_total['alt'] = Mut_total['alt'].str.upper()
Mut_total['mut_info'] = Mut_total['pos'].astype(str) + '_' + Mut_total['ref'] + '_' + Mut_total['alt']

Sample_mut_info_count = Mut_total.value_counts(['SampleName','mut_info']).reset_index()
Sample_mut_info_count.columns = ['sample','mut_info','mut_info_count']

Sample_count = Mut_total.drop_duplicates('SampleName_UMI').value_counts(['SampleName']).reset_index()
Sample_count.columns = ['sample','Sample_count']

Sample_mut_info_count = pd.merge(Sample_mut_info_count,Sample_count,on = 'sample')
Sample_mut_info_count['mut_freq'] = Sample_mut_info_count['mut_info_count'] / Sample_mut_info_count['Sample_count']
Sample_mut_info_count = Sample_mut_info_count.sort_values('mut_info').reset_index(drop = True)
Sample_mut_info_count.to_csv(Mutation_Frequency_file, sep = '\t', index = False)
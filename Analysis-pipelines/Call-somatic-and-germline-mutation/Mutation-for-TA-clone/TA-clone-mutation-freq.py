import pandas as pd

path = '../' # path of the file: Progeny_CallSNP_Plant1.txt

Mut_total = pd.read_table((path + 'Progeny_CallSNP_Plant1.txt'))

Sample_mut_info_count = Mut_total.value_counts(['SampleName','mut_info']).reset_index()
Sample_mut_info_count.columns = ['SampleName','mut_info','mut_info_count']

Sample_count = Mut_total.drop_duplicates('SampleName_UMI').value_counts(['SampleName']).reset_index()
Sample_count.columns = ['SampleName','Sample_count']

Sample_mut_info_count = pd.merge(Sample_mut_info_count,Sample_count,on = 'SampleName')
Sample_mut_info_count['mut_freq'] = Sample_mut_info_count['mut_info_count'] / Sample_mut_info_count['Sample_count']
Sample_mut_info_count = Sample_mut_info_count.sort_values('mut_info').reset_index(drop = True)

Sample_mut_info_count.to_csv((path + 'MutFreq_In_ProgenySample_Plant1.txt'), sep = '\t', index = False)
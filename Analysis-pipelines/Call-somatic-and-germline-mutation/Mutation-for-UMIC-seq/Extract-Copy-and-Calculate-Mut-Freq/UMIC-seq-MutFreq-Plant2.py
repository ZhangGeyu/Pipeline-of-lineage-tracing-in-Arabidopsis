import pandas as pd

HighFreqMut_list_886copy = ['1213_C_+5GTGTG','886_G_A','862_G_C','874_G_A','837_G_C','904_G_T','841_C_T','458_T_-20AACAGGGTAATGAGCCGCAC']

path = '../' # path of the folder contains Parental_CallSNP_Plant2_886copy.txt
copy = '886copy'

Mut_total = pd.read_table((path + 'Parental_CallSNP_Plant2_' + copy + '.txt'))
Mut_total = Mut_total[~Mut_total['mut_info'].isin(HighFreqMut_list_886copy)]

Sample_mut_info_count = Mut_total.value_counts(['SampleName','mut_info']).reset_index()
Sample_mut_info_count.columns = ['SampleName','mut_info','mut_info_count']

Sample_count = Mut_total.drop_duplicates('SampleName_UMI').value_counts(['SampleName']).reset_index()
Sample_count.columns = ['SampleName','Sample_count']

Sample_mut_info_count = pd.merge(Sample_mut_info_count,Sample_count,on = 'SampleName')
Sample_mut_info_count['mut_freq'] = Sample_mut_info_count['mut_info_count'] / Sample_mut_info_count['Sample_count']
Sample_mut_info_count = Sample_mut_info_count.sort_values('mut_info').reset_index(drop = True)

Sample_mut_info_count.to_csv((path + 'MutFreq_In_ParentalSample_Plant2_' + copy + '.txt'), sep = '\t', index = False)
Sample_mut_info_count
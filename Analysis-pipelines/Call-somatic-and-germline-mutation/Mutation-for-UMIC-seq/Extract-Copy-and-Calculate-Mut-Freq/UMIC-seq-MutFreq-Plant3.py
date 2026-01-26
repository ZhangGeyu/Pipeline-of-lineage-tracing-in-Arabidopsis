import pandas as pd

HighFreqMut_list_181copy = ['1213_C_+5GTGCT','181_G_A','182_G_A','183_G_A','184_G_A','736_G_A','601_G_A','657_G_A','719_C_T','238_G_A','658_G_A','663_C_T']

path = '../' # path of the folder contains Parental_CallSNP_Plant3_181copy.txt
copy = '181copy'

Mut_total = pd.read_table((path + 'Parental_CallSNP_Plant3_' + copy + '.txt'))
Mut_total = Mut_total[~Mut_total['mut_info'].isin(HighFreqMut_list_181copy)]

Sample_mut_info_count = Mut_total.value_counts(['SampleName','mut_info']).reset_index()
Sample_mut_info_count.columns = ['SampleName','mut_info','mut_info_count']

Sample_count = Mut_total.drop_duplicates('SampleName_UMI').value_counts(['SampleName']).reset_index()
Sample_count.columns = ['SampleName','Sample_count']

Sample_mut_info_count = pd.merge(Sample_mut_info_count,Sample_count,on = 'SampleName')
Sample_mut_info_count['mut_freq'] = Sample_mut_info_count['mut_info_count'] / Sample_mut_info_count['Sample_count']
Sample_mut_info_count = Sample_mut_info_count.sort_values('mut_info').reset_index(drop = True)

Sample_mut_info_count.to_csv((path + 'MutFreq_In_ParentalSample_Plant3_' + copy + '.txt'), sep = '\t', index = False)
Sample_mut_info_count
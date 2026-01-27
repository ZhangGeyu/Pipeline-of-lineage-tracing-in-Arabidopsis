import pandas as pd

# Pearson's correlation coefficient of rosette leaf mutation frequencies.

Hotspot_Plant1 = ["1092_T_-1C","1114_A_G","1226_C_T","1335_C_T","445_C_T","513_C_G","53_C_-2GA","77_G_A","795_G_T","965_G_A"]

# Raw file
MutFreq_In_ParentalSample = pd.read_table("MutFreq_In_ParentalSample_Plant1.txt")
MutFreq_In_ParentalSample = MutFreq_In_ParentalSample[~MutFreq_In_ParentalSample['mut_info'].isin(Hotspot_Plant1)] # remove hotspot mutations 

# Sample-mutation frequency matrix construction.
Sample_MutFreq_Matrix = MutFreq_In_ParentalSample.pivot_table(index='mut_info', columns='SampleName', values='mut_freq', fill_value=0)

# Pearson's correlation coefficient.
SampleCor = Sample_MutFreq_Matrix.corr(method='pearson').reset_index()

SampleList = list(SampleCor['SampleName'].unique())
sampleList_locals = []
# For each sample, create a correlation dataframe and label groups
for i in SampleList:
    sampleList_locals.append((i+'_corr'))
    locals()[i+'_corr'] = SampleCor[['SampleName',i]]
    locals()[i+'_corr'].columns = ['index','pcc']
    locals()[i+'_corr']['SampleName'] = i
    locals()[i+'_corr']['group'] = ''

    for j in range(len(locals()[i+'_corr']['index'])):
        if locals()[i+'_corr']['index'][j][1] == i[1]:
            locals()[i+'_corr']['group'][j] = 'Within_Branch'
        else:
            locals()[i+'_corr']['group'][j] = 'Between_Branch'

name_1 = sampleList_locals[0]
name_2 = sampleList_locals[1]

CorrMatrix_concat = pd.concat([locals()[name_1],locals()[name_2]])

sampleList_locals.remove(name_1)
sampleList_locals.remove(name_2)

for i in sampleList_locals:
    CorrMatrix_concat = pd.concat([CorrMatrix_concat,locals()[i]])
#  Keep only non-rosette leaf samples (sample name does not start with 'R')
CorrMatrix_concat = CorrMatrix_concat[(CorrMatrix_concat['pcc'] != 1)&\
                                      (CorrMatrix_concat['SampleName'].str.get(0) != 'R')&\
                                        (CorrMatrix_concat['index'].str.get(0) != 'R')].reset_index(drop = True)
CorrMatrix_concat.columns = ['sample1','pcc','sample2','group']

CorrMatrix_concat['combined'] = CorrMatrix_concat[['sample1', 'sample2']].apply(sorted, axis=1).astype(str)
CorrMatrix_concat.drop_duplicates(subset='combined', keep='first', inplace=True)
CorrMatrix_concat.drop('combined', axis=1, inplace=True)
CorrMatrix_concat = CorrMatrix_concat.reset_index(drop = True)

# Clean file
CorrMatrix_concat.to_csv('Sample_MutFreq_correlation.txt'), sep = '\t', index = False)

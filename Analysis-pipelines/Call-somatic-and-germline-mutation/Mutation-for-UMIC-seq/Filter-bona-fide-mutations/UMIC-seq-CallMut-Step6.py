import pandas as pd

path = '../'
depth_cutoff = 10
freq_cutoff = 0.6

Raw_Mutation = pd.read_table(path + 'ConsensusSequence_CallSNP_Raw.txt')

Raw_Mutation = Raw_Mutation[Raw_Mutation['depth'] >= depth_cutoff].reset_index(drop = True)
Raw_Mutation = Raw_Mutation[Raw_Mutation['freq'] >= freq_cutoff].reset_index(drop = True)
Raw_Mutation['UMI'] = Raw_Mutation['Sample_UMI'].str.split('_').str.get(1)
Raw_Mutation['SampleName_UMI'] = Raw_Mutation['SampleName'] + '_' + Raw_Mutation['UMI']
Raw_Mutation['mut_info'] = Raw_Mutation['pos'].astype(str) + '_' + Raw_Mutation['ref'] + '_' + Raw_Mutation['alt']

Raw_Mutation = Raw_Mutation[['pos','ref','alt','mut_info','SampleName','UMI','SampleName_UMI']]

Raw_Mutation.to_csv((path + 'ConsensusSequence_CallSNP.txt'), sep = '\t')
import pandas as pd

Folder_contain_mpileup_file = '/path/to/your/folder' # same as Output_Folder_for_Splited_Bam in 5.Call_mutations_from_Readouts
SampleBC_SampleName_file = 'SampleBC_SampleName.txt'
Raw_Mutation_file = ''
Clean_Mutation_file = ''

depth_cutoff = 10
freq_cutoff = 0.6

BC_name = pd.read_table(SampleBC_SampleName_file)
BC_name_list = list(BC_name['SampleBC'].unique())

mpileup_files = [f.split('.')[0] for f in os.listdir(Folder_contain_mpileup_file) if f.endswith('.mpileup')]
mpileup_files = [f for f in mpileup_files if f.split('_')[0] in BC_name_list]

alt_count_df = pd.DataFrame(columns=['pos','ref','alt','ReadsName'])

def replace_with_pipe(match):
    before, n, after = match.groups()
    n = int(n)
    return before + str(n) + after[:n] + '|' + after[n:]

## Extract all mutations in the .mpileup files

with open(Raw_Mutation_file,'w') as alt_count_df:
    alt_count_df.write('pos\t'+'ref\t'+'alt\t'+'Sample_UMI\t'+'SampleName\t'+'alt_count\t'+'depth\t'+'freq\n')

    n = 0
    for m_file in mpileup_files:
        n = n + 1
        if n % 1000 == 0:
            print(m_file.split('_')[0] + ': ' + str(n) + ' / ' + str(len(mpileup_files)))
        SampleName = BC_name[BC_name['SampleBC'] == m_file.split('_')[0]].reset_index()['SampleName'][0]
        mpileup_split = pd.read_table((Folder_contain_mpileup_file + m_file + '.mpileup'),\
                                    names = ['chr','pos', 'ref','depth','alt','quality','?']).reset_index()
        mpileup_split = mpileup_split[['pos','ref','depth','alt']].reset_index(drop = True)

        for i in range(len(mpileup_split['pos'])):
            pos = mpileup_split['pos'][i]
            ref = mpileup_split['ref'][i].upper()
            alt = mpileup_split['alt'][i].upper()
            pattern_1 = r'(\D)(\d+)([acgtnACGTN]+)'
            alt = re.sub(pattern_1, replace_with_pipe, alt)
            alt = re.sub('\^.','',alt)
            alt = re.sub('\$','',alt)
            pattern = r'([,.][+-]\d+[ACGTNacgtn]+)|([+-]\d+[ACGTNacgtn]+)|([,.<>*])|([ACGTNacgtn])'
            matches = re.findall(pattern, alt)
            alt_list = [item for tup in matches for item in tup if item != '']
            alt_list = [item[1:] if len(item) > 1 and (item[0] == '.' or item[0] == ',') else item for item in alt_list]
            alt_list = [s.replace(',', '.').replace('*', '.').replace('>', '.').replace('<', '.') for s in alt_list]
            #print(alt_list)
            element_counts = dict(Counter(alt_list))
            #print(element_counts)
            for alt_i in element_counts:
                if alt_i != '.':
                    freq = element_counts[alt_i] / mpileup_split['depth'][i]
                    alt_count_df.write(pos.astype(str) + '\t' + ref + '\t'+ alt_i + '\t' + m_file + '\t' + SampleName + '\t' + \
                        str(element_counts[alt_i]) + '\t' + str(mpileup_split['depth'][i]) + '\t' + str(freq) + '\n')

## Set the 0.6 cutoff for all the mutations to distinguish the true mutations on each readout

Raw_Mutation = pd.read_table(Raw_Mutation_file)
Raw_Mutation = Raw_Mutation[Raw_Mutation['depth'] >= depth_cutoff].reset_index(drop = True)
All_Sample_UMI_df = Raw_Mutation[['Sample_UMI','SampleName']].drop_duplicates('Sample_UMI').reset_index(drop = True)
Raw_Mutation = Raw_Mutation[Raw_Mutation['freq'] >= freq_cutoff].reset_index(drop = True)

Raw_Mutation['UMI'] = Raw_Mutation['Sample_UMI'].str.split('_').str.get(1)
Raw_Mutation['SampleName_UMI'] = Raw_Mutation['SampleName'] + '_' + Raw_Mutation['UMI']
Raw_Mutation.to_csv(Clean_Mutation_file, sep = '\t')

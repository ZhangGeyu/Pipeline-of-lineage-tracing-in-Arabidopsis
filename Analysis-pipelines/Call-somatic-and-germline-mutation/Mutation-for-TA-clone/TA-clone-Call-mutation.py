import re
import pandas as pd

path = '../' # Path of the Plant1_Progeny_merge_seq.mpileup

mpileup_file = pd.read_table((path + 'Plant1_Progeny_merge_seq.mpileup'),\
    names = ['pos', 'ref','depth','alt','quality','?','ReadsName'])
mpileup_file = mpileup_file[mpileup_file['depth'] != '0'].reset_index().drop('?', axis = 1)
mpileup_file.columns = ['chr','pos','ref','depth','alt','quality','ReadsName']

mpileup_file_HaveMut = mpileup_file[mpileup_file['alt'].str.contains('[a-zA-Z]')].reset_index(drop = True)

def replace_with_pipe(match):
    before, n, after = match.groups()
    n = int(n)
    return before + str(n) + after[:n] + '|' + after[n:]

Mut_reads_unique = pd.DataFrame(columns=['chr', 'pos', 'ref','alt','ReadsName'])

for i in range(len(mpileup_file_HaveMut['chr'])):
    #print(i)
    ref_chr = mpileup_file_HaveMut['chr'][i]
    ref = mpileup_file_HaveMut['ref'][i].upper()
    pos = mpileup_file_HaveMut['pos'][i]
    alt = mpileup_file_HaveMut['alt'][i].upper()
    
    pattern_1 = r'(\D)(\d+)([acgtnACGTN]+)'
    alt = re.sub(pattern_1, replace_with_pipe, alt)
    alt = re.sub('\^.','',alt)
    alt = re.sub('\$','',alt)
    pattern = r'([,.][+-]\d+[ACGTNacgtn]+)|([+-]\d+[ACGTNacgtn]+)|([,.<>*])|([ACGTNacgtn])'
    matches = re.findall(pattern, alt)
    alt_list = [item for tup in matches for item in tup if item != '']
    alt_list = [item[1:] if len(item) > 1 and (item[0] == '.' or item[0] == ',') else item for item in alt_list]
    alt_list = [s.replace(',', '.').replace('*', '.').replace('>', '.').replace('<', '.') for s in alt_list]

    ReadsName_list = mpileup_file_HaveMut['ReadsName'][i].split(',')

    if len(ReadsName_list) == len(alt_list):
        for j in range(len(ReadsName_list)):
            if alt_list[j] != '.':
                Mut_reads_unique = Mut_reads_unique.append({'chr':ref_chr,'pos':pos,'ref': ref,\
                'alt': alt_list[j],'ReadsName':ReadsName_list[j]}, ignore_index=True)
            else:
                continue

Mut_reads_unique['SampleName'] = Mut_reads_unique['ReadsName'].str.rsplit('-', n=1).str[0]
Mut_reads_unique['clone'] = Mut_reads_unique['ReadsName'].str.split('-').str.get(-1)
Mut_reads_unique['SampleName_UMI'] = Mut_reads_unique['SampleName'] + '_' + Mut_reads_unique['clone']
Mut_reads_unique['mut_info'] = Mut_reads_unique['pos'].astype(str) + '_' + Mut_reads_unique['ref'] + '_' + Mut_reads_unique['alt']

Mut_reads_unique.to_csv((path + 'Progeny_CallSNP_Plant1.txt'), sep = '\t', index = False)
Mut_reads_unique
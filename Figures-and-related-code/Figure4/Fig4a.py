# Extract parental and germline sequences to create a FASTA file for phylogenetic tree construction.

from Bio import SeqIO
import pandas as pd
import random
import csv

Hotspot_Plant2 = [""]
# remove copy mutations
HighFreqMut_list_886copy = ['1213_C_+5GTGTG','886_G_A','862_G_C','874_G_A','837_G_C','904_G_T','841_C_T','458_T_-20AACAGGGTAATGAGCCGCAC'] 

# Offspring

Progeny_mutation = pd.read_table('MutFreq_In_ProgenySample_Plant2_886copy.txt')
Progeny_mutation = Progeny_mutation[(Progeny_mutation['mut_freq'] >= 0.5)&\
    (~Progeny_mutation['mut_info'].isin(Hotspot_Plant2))&(~Progeny_mutation['mut_info'].isin(HighFreqMut_list_886copy))] 
#  Extract position, reference and alternative bases from mutation info
Progeny_mutation['pos'] = Progeny_mutation['mut_info'].str.split('_').str.get(0).astype(int)
Progeny_mutation['ref'] = Progeny_mutation['mut_info'].str.split('_').str.get(1)
Progeny_mutation['alt'] = Progeny_mutation['mut_info'].str.split('_').str.get(2)
Progeny_mutation['SampleName_UMI'] = Progeny_mutation['SampleName'] + '_progeny'

# Parental
# Read parental mutation data
Parental_mutation = pd.read_table('Parental_CallSNP_Plant2_886copy.txt')
Parental_mutation = Parental_mutation[(~Parental_mutation['mut_info'].isin(Hotspot_Plant2))&(~Parental_mutation['mut_info'].isin(HighFreqMut_list_886copy))]

print('Number of total Readout: ' + str(Parental_mutation['SampleName_UMI'].nunique()))
# Count mutations per UMI
Sample_UMI_MutCount = Parental_mutation.value_counts('SampleName_UMI').reset_index()
Sample_UMI_MutCount.columns = ['SampleName_UMI','count']

UMI_MutCount_list = list(Sample_UMI_MutCount[Sample_UMI_MutCount['count'] <= Sample_UMI_MutCount['count'].mean() + Sample_UMI_MutCount['count'].std()]['SampleName_UMI'].unique())
Parental_mutation = Parental_mutation[Parental_mutation['SampleName_UMI'].isin(UMI_MutCount_list)] # UMI MutCount cutoff


print('Number of Readout with suitable mut count: ' + str(Parental_mutation['SampleName_UMI'].nunique()))
# For each parental sample, randomly select up to 50 UMIs to reduce data size
for SampleName in Parental_mutation['SampleName'].unique():
    SampleName_UMI_list = list(Parental_mutation[Parental_mutation['SampleName'] == SampleName]['SampleName_UMI'].unique())
    if len(SampleName_UMI_list) >= 50:
        Not_SampleName_UMI_list = set(SampleName_UMI_list) - set(random.sample(SampleName_UMI_list,50)) 
        Parental_mutation = Parental_mutation[~Parental_mutation['SampleName_UMI'].isin(Not_SampleName_UMI_list)]
    else:
        continue

print('Parental random select Readout count: ' + str(Parental_mutation['SampleName_UMI'].nunique()))
#  Merge parental and progeny mutation data
UMI_mut_total = pd.concat([Progeny_mutation,Parental_mutation]).reset_index(drop = True)
UMI_mut_total[['SampleName_UMI','mut_info']].to_csv('ParentalAndOffspring_mut.txt', sep = '\t', index = False)

# Generate consensus sequence

reference = SeqIO.read('reference-Plant2.fa', "fasta")

Cutoff_Mutation = UMI_mut_total.copy()
BC_name_list = Cutoff_Mutation['SampleName_UMI'].unique()

rows = []
# Generate consensus sequence for each UMI
for BC in BC_name_list:
    refstr = list(str(reference.seq).upper())
    mpileup_split_alt = Cutoff_Mutation[Cutoff_Mutation['SampleName_UMI'] == BC].reset_index(drop = True)

    for i in range(len(mpileup_split_alt['pos'])):
        pos = mpileup_split_alt['pos'][i] - 1
        ref = mpileup_split_alt['ref'][i]
        alt = mpileup_split_alt['alt'][i]
        SampleName = mpileup_split_alt['SampleName'][i]
        if refstr[pos] == ref:
            if '+' in alt:
                refstr[pos] = refstr[pos] + ''.join(char for char in alt if char.isalpha())
            elif '-' in alt:
                digits = int(''.join(char for char in alt if char.isdigit()))
                for j in range(1,digits + 1):
                    refstr[pos+j] = ''
            else:
                refstr[pos] = alt

    consensus_seq = ''.join(refstr)
    consensus_seq = consensus_seq.upper()
    UMI = BC.split('_')[1]
    rows.append({'SampleName':SampleName, 'seq':consensus_seq, 'UMI':UMI, 'SampleName_UMI':BC})

ConsensueSeq = pd.DataFrame(rows)
# Write consensus sequences to FASTA file
with open('ConsensusSeq_Plant2.fasta', 'w') as fasta_consensus:
    for index, row in ConsensueSeq.iterrows():
        reads_id = row['SampleName_UMI'].replace('-','_')
        seq = row['seq']
        fasta_consensus.write(f'>{reads_id}\n{seq}\n')
    fasta_consensus.write(f'>reference\n{str(reference.seq).upper()}\n')
# Count UMIs in different branches
data = list(ConsensueSeq['SampleName_UMI'].unique())
count_1 = sum(s.startswith('B1') for s in data)
count_2 = sum(s.startswith('B2') for s in data)
count_3 = sum(s.startswith('B3') for s in data)
count_R = sum(s.startswith('R') for s in data)

print("Branch1 UMI count: ", count_1)
print("Branch2 UMI count: ", count_2)
print("Branch3 UMI count: ", count_3)
print("Rosette leaves UMI count: ", count_R)

# Tree colour annotation

header = ['DATASET_COLORSTRIP','SEPARATOR TAB','DATASET_LABEL	Region1 Sample Colors','COLOR	#808080',\
'LEGEND_TITLE	Sample','LEGEND_SHAPES	1	1	1	1	1	1','LEGEND_COLORS	#1ec0ff	#9381ff	#3ac569	#fb8500	#D3D3D3	#D3D3D3',\
'LEGEND_LABELS	5d	8d	10d	14d	16d	18d','STRIP_WIDTH	25','MARGIN	5','SHOW_INTERNAL	0','DATA']

UMI_mut_total = pd.read_table('ParentalAndOffspring_mut.txt')
UMI_mut_total = UMI_mut_total.drop_duplicates('SampleName_UMI').reset_index(drop = True)
UMI_mut_total['SampleName_UMI'] = UMI_mut_total['SampleName_UMI'].str.replace('-','_')

UMI_mut_total.loc[UMI_mut_total['SampleName_UMI'].str.contains('P'), 'group'] = 'Offspring_Haplotype'
UMI_mut_total.loc[UMI_mut_total['SampleName_UMI'].str.contains('CL'), 'group'] = 'Parental_Leaves'
UMI_mut_total.loc[UMI_mut_total['SampleName_UMI'].str.contains('RL'), 'group'] = 'Parental_Rosette_Leaves'

UMI_mut_total['Branch'] = UMI_mut_total['SampleName_UMI'].str.get(1)

UMI_mut_total['type'] = 'range'
UMI_mut_total['colour'] = '#D3D3D3'

UMI_mut_total.loc[(UMI_mut_total['group'] == 'Offspring_Haplotype')&(UMI_mut_total['Branch'] == '1'), 'colour'] = '#1EC0FF'
UMI_mut_total.loc[(UMI_mut_total['group'] == 'Offspring_Haplotype')&(UMI_mut_total['Branch'] == '2'), 'colour'] = '#9381FF'
UMI_mut_total.loc[(UMI_mut_total['group'] == 'Offspring_Haplotype')&(UMI_mut_total['Branch'] == '3'), 'colour'] = '#3AC569'
UMI_mut_total.loc[(UMI_mut_total['group'] == 'Offspring_Haplotype')&(UMI_mut_total['Branch'] == '4'), 'colour'] = '#ee9b00'

UMI_mut_total['normal'] = 'normal'
UMI_mut_total['width'] = 2

with open('Tree_annotation_G.txt', mode='w', newline='') as file:
    writer = csv.writer(file)
    for line in header:
        writer.writerow([line])
    UMI_mut_total[['SampleName_UMI','colour']].to_csv(file, sep = '\t',index = False,header = False)

UMI_mut_total['type'] = 'range'
UMI_mut_total['colour'] = '#D3D3D3'
UMI_mut_total.loc[(UMI_mut_total['group'] == 'Parental_Leaves')&(UMI_mut_total['Branch'] == '1'), 'colour'] = '#1EC0FF'
UMI_mut_total.loc[(UMI_mut_total['group'] == 'Parental_Leaves')&(UMI_mut_total['Branch'] == '2'), 'colour'] = '#9381FF'
UMI_mut_total.loc[(UMI_mut_total['group'] == 'Parental_Leaves')&(UMI_mut_total['Branch'] == '3'), 'colour'] = '#3AC569'
UMI_mut_total.loc[(UMI_mut_total['group'] == 'Parental_Leaves')&(UMI_mut_total['Branch'] == '4'), 'colour'] = '#ee9b00'
UMI_mut_total.loc[(UMI_mut_total['group'] == 'Parental_Rosette_Leaves'), 'colour'] = '#585858'

UMI_mut_total['normal'] = 'normal'
UMI_mut_total['width'] = 2

with open('Tree_annotation_P.txt', mode='w', newline='') as file:
    writer = csv.writer(file)
    for line in header:
        writer.writerow([line])
    UMI_mut_total[['SampleName_UMI','colour']].to_csv(file, sep = '\t',index = False,header = False)

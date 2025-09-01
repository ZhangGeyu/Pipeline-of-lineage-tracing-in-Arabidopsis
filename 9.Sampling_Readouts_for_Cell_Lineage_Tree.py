### 子代 UMI + TA clone + 亲代 seeds + 亲代其它

# Plant 2
HighFreqMut_list = ['886_G_A','862_G_C','874_G_A','837_G_C','904_G_T','841_C_T','458_T_-20AACAGGGTAATGAGCCGCAC']
# Plant 3
HighFreqMut_list = ['181_G_A','182_G_A','183_G_A','184_G_A','736_G_A','601_G_A','657_G_A','719_C_T','238_G_A','658_G_A','663_C_T']


Offspring_Haplotype_cutoff = 0.5
Mutation_file_For_Copy_Used_Progeny = ''
Mutation_file_For_Copy_Used_Somatic = ''
reference_file = '/path/to/your/folder/reference.fa'
Sequence_txt_file = ''
Sequence_fasta_file = ''

Offspring_UMI_mut = pd.read_table(Mutation_file_For_Copy_Used_Progeny)
Offspring_UMI_mut = Offspring_UMI_mut[Offspring_UMI_mut['mut_freq'] >= Offspring_Haplotype_cutoff]

Offspring_UMI_mut['SampleName_UMI'] = Offspring_UMI_mut['sample']
Offspring_UMI_mut['pos'] = Offspring_UMI_mut['mut_info'].str.split('_').str.get(0).astype(int)
Offspring_UMI_mut['ref'] = Offspring_UMI_mut['mut_info'].str.split('_').str.get(1)
Offspring_UMI_mut['alt'] = Offspring_UMI_mut['mut_info'].str.split('_').str.get(2)
Offspring_UMI_mut = Offspring_UMI_mut[['pos','ref','alt','mut_info','sample','UMI','SampleName_UMI']]
Offspring_UMI_mut.columns = ['pos','ref','alt','mut_info','SampleName','UMI','SampleName_UMI']
Offspring_UMI_mut['group'] = 'Offspring_Haplotype'
Offspring_UMI_mut['Branch'] = Offspring_UMI_mut['SampleName'].str.split('-').str.get(1)
Offspring_UMI_mut = Offspring_UMI_mut[(Offspring_UMI_mut['pos'] != 1213)]

Sample_UMI_MutCount = Offspring_UMI_mut.value_counts('SampleName_UMI').reset_index()
Sample_UMI_MutCount.columns = ['SampleName_UMI','count']

UMI_MutCount_cutoff = len(HighFreqMut_list) + 1
UMI_MutCount_list = list(Sample_UMI_MutCount[Sample_UMI_MutCount['count'] >= UMI_MutCount_cutoff]['SampleName_UMI'].unique())
Offspring_UMI_mut = Offspring_UMI_mut[Offspring_UMI_mut['SampleName_UMI'].isin(UMI_MutCount_list)] # UMI MutCount cutoff


# Parental 
Parental_UMI_mut = pd.read_table(Mutation_file_For_Copy_Used_Somatic, usecols = ['pos','ref','alt','SampleName','UMI','SampleName_UMI'])
Parental_UMI_mut['mut_info'] = Parental_UMI_mut['pos'].astype(str) + '_' + Parental_UMI_mut['ref'] + '_' + Parental_UMI_mut['alt']
Parental_UMI_mut = Parental_UMI_mut[(Parental_UMI_mut['pos'] != 1213)]

Parental_UMI_mut['Branch'] = Parental_UMI_mut['SampleName'].str.split('-').str.get(0)
Parental_UMI_mut.loc[Parental_UMI_mut['SampleName_UMI'].str.contains('RL'), 'group'] = 'Parental_Rosette_Leaves'
Parental_UMI_mut.loc[Parental_UMI_mut['SampleName_UMI'].str.contains('CL'), 'group'] = 'Parental_Leaves'

Sample_UMI_MutCount = Parental_UMI_mut.value_counts('SampleName_UMI').reset_index()
Sample_UMI_MutCount.columns = ['SampleName_UMI','count']
UMI_MutCount_cutoff = len(HighFreqMut_list) + 1
UMI_MutCount_list = list(Sample_UMI_MutCount[Sample_UMI_MutCount['count'] >= UMI_MutCount_cutoff]['SampleName_UMI'].unique())
Parental_UMI_mut = Parental_UMI_mut[Parental_UMI_mut['SampleName_UMI'].isin(UMI_MutCount_list)] # UMI MutCount cutoff

for SampleName in Parental_UMI_mut['SampleName'].unique():
    SampleName_UMI_list = list(Parental_UMI_mut[Parental_UMI_mut['SampleName'] == SampleName]['SampleName_UMI'].unique())
    if len(SampleName_UMI_list) >= 50:
        Not_SampleName_UMI_list = set(SampleName_UMI_list) - set(random.sample(SampleName_UMI_list,50))
        Parental_UMI_mut = Parental_UMI_mut[~Parental_UMI_mut['SampleName_UMI'].isin(Not_SampleName_UMI_list)]
    else:
        continue

UMI_mut_total = pd.concat([Offspring_UMI_mut,Parental_UMI_mut]).reset_index(drop = True)

# Generate consensus sequence

reference = SeqIO.read(reference_file, "fasta")

Cutoff_Mutation = UMI_mut_total.copy()

BC_name_list = Cutoff_Mutation['SampleName_UMI'].unique()
print(len(BC_name_list))

with open(Sequence_txt_file,'w') as ConsensueSeq:
    ConsensueSeq.write('SampleName\t'+'seq\t'+'UMI\t'+'SampleName_UMI\n')

    n = 0
    for BC in BC_name_list:
        n = n + 1
        if n % 1000 == 0:
            print('Plant:9 ' + str(n) + ': ' + BC)
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
        ConsensueSeq.write(SampleName + '\t'+ consensus_seq + '\t' + UMI + '\t' + BC + '\n')

ConsensueSeq = pd.read_table(Sequence_txt_file, usecols = ['SampleName_UMI','seq'])
ConsensueSeq['SampleName_UMI'] = ConsensueSeq['SampleName_UMI'].str.replace('(','.').str.replace(')','.')

with open(Sequence_fasta_file, 'w') as fasta_consensus:
    for index, row in ConsensueSeq.iterrows():
        reads_id = row['SampleName_UMI'].replace('-','_')
        seq = row['seq']
        fasta_consensus.write(f'>{reads_id}\n{seq}\n')
    fasta_consensus.write(f'>reference\n{str(reference.seq).upper()}\n')

## Command for Mutiple Alignment
'''
nohup mafft --thread 50 ConsensusSequence_181copy.fasta > ConsensusSequence_181copy_Multiple_Alignment.fasta &
'''
## Command for Tree Construction
'''
nohup FastTree -nt ConsensusSequence_181copy_Multiple_Alignment.fasta > ConsensusSequence_181copy_Tree.nwk &
'''

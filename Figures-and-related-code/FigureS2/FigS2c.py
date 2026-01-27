Hotspot_Plant1 = ["1092_T_-1C","1114_A_G","1226_C_T","1335_C_T","445_C_T","513_C_G","53_C_-2GA","77_G_A","795_G_T","965_G_A"]

# Raw file
Parental_mutation = pd.read_table('Parental_CallSNP_Plant1.txt')
Parental_mutation = Parental_mutation[~Parental_mutation['mut_info'].isin(Hotspot_Plant1)]
print(len(Parental_mutation['SampleName_UMI'].unique()))

## Generate Consensus sequence
reference = SeqIO.read('reference-Plant1.fa', "fasta")

BC_name_list = Parental_mutation['SampleName_UMI'].unique()   # Get list of unique UMI names
rows = []

n = 0
for BC in BC_name_list:
    n = n + 1
    if n % 1000 == 0:
        print('Plant1: ' + str(n) + ': ' + BC)
    refstr = list(str(reference.seq).upper())
    mpileup_split_alt = Parental_mutation[Parental_mutation['SampleName_UMI'] == BC].reset_index(drop = True)    #Get all mutations for current UMI sample

    for i in range(len(mpileup_split_alt['pos'])):
        pos = mpileup_split_alt['pos'][i] - 1
        ref = mpileup_split_alt['ref'][i]
        alt = mpileup_split_alt['alt'][i]
        SampleName = mpileup_split_alt['SampleName'][i]
        if refstr[pos] == ref:
            if '+' in alt:      # Handle insertion mutations
                refstr[pos] = refstr[pos] + ''.join(char for char in alt if char.isalpha())
            elif '-' in alt:    # Handle deletion mutations
                digits = int(''.join(char for char in alt if char.isdigit()))
                for j in range(1,digits + 1):
                    refstr[pos+j] = ''
            else:
                refstr[pos] = alt

    consensus_seq = ''.join(refstr).upper()
    UMI = BC.split('_')[1]
    rows.append({'SampleName': SampleName,'seq': consensus_seq,'UMI': UMI,'SampleName_UMI': BC})

ConsensusSeq_df = pd.DataFrame(rows)    # Convert rows list to DataFrame

### Fraction of # of UMI In each Allele
# Count occurrences of each sequence
Allele_UMI_Count = ConsensusSeq_df.value_counts(['seq']).reset_index()
Allele_UMI_Count.columns = ['seq','count']

Allele_UMI_Count['Redundancy'] = ''
Allele_UMI_Count.loc[Allele_UMI_Count['count'] < 5, 'Redundancy'] = Allele_UMI_Count['count']
Allele_UMI_Count.loc[Allele_UMI_Count['count'] >= 5, 'Redundancy'] = '>=5'    # For counts >= 5, mark as >=5

Allele_Redundancy_Count = Allele_UMI_Count.value_counts('Redundancy').reset_index()
Allele_Redundancy_Count.columns = ['Redundancy','count']
Allele_Redundancy_Count['fraction'] = Allele_Redundancy_Count['count'] / sum(Allele_Redundancy_Count['count'])   #Calculate fraction for each category

# Clean file
Allele_Redundancy_Count.to_csv('Allele_Redundancy_Count.txt', sep = '\t')

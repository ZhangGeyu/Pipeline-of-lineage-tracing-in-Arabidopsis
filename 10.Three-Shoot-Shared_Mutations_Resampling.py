Offspring_Haplotype_cutoff = 0.5
Mutation_Frequency_file_Progeny = ''
Mutation_file_For_Copy_Used_Somatic = ''
Branch_Shared_Mutation_file = ''
Readout_Mutation_Count_file = ''


Offspring_Mut = pd.read_table(Mutation_Frequency_file_Progeny), usecols = ['sample','mut_info','mut_freq'])
Offspring_Mut.columns = ['SampleName','mut_info','freq']
Offspring_Mut = Offspring_Mut[(Offspring_Mut['freq'] >= 0.5)]

Offspring_Mut['Branch'] = Offspring_Mut['SampleName'].str.get(1)

Parental_Mut = pd.read_table(Mutation_file_For_Copy_Used_Somatic)
Parental_Mut['mut_info'] = Parental_Mut['pos'].astype(str) + '_' + Parental_Mut['ref'].str.upper() + '_' + Parental_Mut['alt'].str.upper()


Parental_Mut['UMI'] = Parental_Mut['SampleName_UMI'].str.split('_').str.get(-1)
Parental_Mut['SampleName_UMI'] = Parental_Mut['SampleName'] + '_' + Parental_Mut['UMI']
Parental_Mut['Branch'] = Parental_Mut['SampleName'].str.get(1)

Branch_list = list(set(Parental_Mut['Branch'].unique()) & set(Offspring_Mut['Branch'].unique()))
print(Branch_list)

Offspring_Mut = Offspring_Mut[Offspring_Mut['Branch'].isin(Branch_list)][['SampleName','mut_info','Branch']]
Parental_Mut = Parental_Mut[Parental_Mut['Branch'].isin(Branch_list)][['SampleName','SampleName_UMI','mut_info','Branch']]

for Branch in Branch_list:
    locals()['Parental_Mut' + Branch] = Parental_Mut[Parental_Mut['Branch'] == Branch]
    locals()['Offspring_Mut' + Branch] = Offspring_Mut[Offspring_Mut['Branch'] == Branch]
    OP_Mut_list = set(locals()['Offspring_Mut' + Branch]['mut_info'].unique())
    locals()['Parental_Mut' + Branch] = locals()['Parental_Mut' + Branch][locals()['Parental_Mut' + Branch]['mut_info'].isin(OP_Mut_list)]

Parental_Mut = pd.concat([locals()['Parental_Mut' + Branch_list[0]],locals()['Parental_Mut' + Branch_list[1]],\
    locals()['Parental_Mut' + Branch_list[2]]]).reset_index(drop = True)
Offspring_Mut = pd.concat([locals()['Offspring_Mut' + Branch_list[0]],locals()['Offspring_Mut' + Branch_list[1]],\
    locals()['Offspring_Mut' + Branch_list[2]]]).reset_index(drop = True)

Offspring_SampleMutCount = Offspring_Mut.value_counts('SampleName').reset_index(name = 'mut_count')
Offspring_SampleMutCount['Branch'] = Offspring_SampleMutCount['SampleName'].str.get(1)

Offspring_SampleMutCount['Group'] = 'Offspring'
Offspring_SampleMutCount.loc[(Offspring_SampleMutCount['mut_count'] >= 1) & (Offspring_SampleMutCount['mut_count'] <= 10), 'Range'] = 'Range1'
Offspring_SampleMutCount.loc[(Offspring_SampleMutCount['mut_count'] >= 11) & (Offspring_SampleMutCount['mut_count'] <= 20), 'Range'] = 'Range2'
Offspring_SampleMutCount.loc[(Offspring_SampleMutCount['mut_count'] >= 21) & (Offspring_SampleMutCount['mut_count'] <= 30), 'Range'] = 'Range3'
Offspring_SampleMutCount.loc[(Offspring_SampleMutCount['mut_count'] >= 31) & (Offspring_SampleMutCount['mut_count'] <= 40), 'Range'] = 'Range4'
Offspring_SampleMutCount.loc[(Offspring_SampleMutCount['mut_count'] >= 41) & (Offspring_SampleMutCount['mut_count'] <= 50), 'Range'] = 'Range5'
Offspring_SampleMutCount.loc[(Offspring_SampleMutCount['mut_count'] >= 51), 'Range'] = 'Range6'

Parental_SampleMutCount = Parental_Mut.value_counts('SampleName_UMI').reset_index(name = 'mut_count')
Parental_SampleMutCount.columns = ['SampleName','mut_count']
Parental_SampleMutCount['Branch'] = Parental_SampleMutCount['SampleName'].str.get(1)

Parental_SampleMutCount['Group'] = 'Parental'
Parental_SampleMutCount.loc[(Parental_SampleMutCount['mut_count'] >= 1) & (Parental_SampleMutCount['mut_count'] <= 10), 'Range'] = 'Range1'
Parental_SampleMutCount.loc[(Parental_SampleMutCount['mut_count'] >= 11) & (Parental_SampleMutCount['mut_count'] <= 20), 'Range'] = 'Range2'
Parental_SampleMutCount.loc[(Parental_SampleMutCount['mut_count'] >= 21) & (Parental_SampleMutCount['mut_count'] <= 30), 'Range'] = 'Range3'
Parental_SampleMutCount.loc[(Parental_SampleMutCount['mut_count'] >= 31) & (Parental_SampleMutCount['mut_count'] <= 40), 'Range'] = 'Range4'
Parental_SampleMutCount.loc[(Parental_SampleMutCount['mut_count'] >= 41) & (Parental_SampleMutCount['mut_count'] <= 50), 'Range'] = 'Range5'
Parental_SampleMutCount.loc[(Parental_SampleMutCount['mut_count'] >= 51), 'Range'] = 'Range6'

SampleMutCount = pd.concat([Offspring_SampleMutCount,Parental_SampleMutCount])

with open(Branch_Shared_Mutation_file,'w') as MutRandom:
    MutRandom.write('Branch\t'+'Group\t'+'Mut_InDiffBranch_Count\n')

    with open(Readout_Mutation_Count_file,'w') as UMI_MutCount:
        UMI_MutCount.write('mut_count\t'+'Group\t'+'Branch\t'+'num\n')
            
        Offspring_Mut_Branch_1 = Offspring_Mut[Offspring_Mut['Branch'] == '1']
        Offspring_Mut_Branch_2 = Offspring_Mut[Offspring_Mut['Branch'] == '2']
        Offspring_Mut_Branch_3 = Offspring_Mut[Offspring_Mut['Branch'] == '3']
        MutOffspring_Mut_BranchCount = len(set(Offspring_Mut_Branch_1['mut_info'].unique()) & set(Offspring_Mut_Branch_2['mut_info'].unique()) & set(Offspring_Mut_Branch_3['mut_info'].unique()))

        Offspring_UMI_MutCount_1 = Offspring_Mut_Branch_1.value_counts('SampleName').reset_index(name = 'mut_count')[['mut_count']]
        Offspring_UMI_MutCount_1['group'] = 'Offspring'
        Offspring_UMI_MutCount_1['Branch'] = 'Branch_1'
        Offspring_UMI_MutCount_1['num'] = 1
        Offspring_UMI_MutCount_1.to_csv(UMI_MutCount, sep = '\t', index=False, header=False)
        Offspring_UMI_MutCount_2 = Offspring_Mut_Branch_2.value_counts('SampleName').reset_index(name = 'mut_count')[['mut_count']]
        Offspring_UMI_MutCount_2['group'] = 'Offspring'
        Offspring_UMI_MutCount_2['Branch'] = 'Branch_2'
        Offspring_UMI_MutCount_2['num'] = 1
        Offspring_UMI_MutCount_2.to_csv(UMI_MutCount, sep = '\t', index=False, header=False)
        Offspring_UMI_MutCount_3 = Offspring_Mut_Branch_3.value_counts('SampleName').reset_index(name = 'mut_count')[['mut_count']]
        Offspring_UMI_MutCount_3['group'] = 'Offspring'
        Offspring_UMI_MutCount_3['Branch'] = 'Branch_3'
        Offspring_UMI_MutCount_3['num'] = 1
        Offspring_UMI_MutCount_3.to_csv(UMI_MutCount, sep = '\t', index=False, header=False)

        n = 1000
        for i in range(1,n+1):
            if i % 500 == 0:
                print(i)
            for Branch in ['1','2','3']:
                Offspring_Mut_subset = Offspring_Mut[Offspring_Mut['Branch'] == Branch]
                Parental_Mut_subset = Parental_Mut[Parental_Mut['Branch'] == Branch]

                Offspring_SampleMutCount_Branch = Offspring_SampleMutCount[Offspring_SampleMutCount['Branch'] == Branch]
                Parental_SampleMutCount_Branch = Parental_SampleMutCount[Parental_SampleMutCount['Branch'] == Branch]

                UMI_Select = []

                for Ra in Offspring_SampleMutCount_Branch['Range'].unique():
                    Offspring_SampleMutCount_RangeCount = len(Offspring_SampleMutCount_Branch[Offspring_SampleMutCount_Branch['Range'] == Ra]['SampleName'])
                    Parental_SampleMutCount_Branch = Parental_SampleMutCount_Branch[~Parental_SampleMutCount_Branch['SampleName'].isin(UMI_Select)]
                    Parental_Range_UMI = list(Parental_SampleMutCount_Branch[Parental_SampleMutCount_Branch['Range'] == Ra]['SampleName'].unique())
                    Ra_Add = Ra
                    
                    while len(Parental_Range_UMI) < Offspring_SampleMutCount_RangeCount:
                        Parental_Range_UMI_add = list(Parental_SampleMutCount_Branch[Parental_SampleMutCount_Branch['Range'] == Ra_Add[:-1] + str(int(Ra_Add[-1]) - 1)]['SampleName'].unique())
                        Parental_Range_UMI = list(set(Parental_Range_UMI) | set(Parental_Range_UMI_add))
                        Ra_Add = Ra_Add[:-1] + str(int(Ra_Add[-1]) - 1)

                        if len(Parental_Range_UMI) >= Offspring_SampleMutCount_RangeCount:
                            break
                    
                    locals()['UMI_' + Ra] = random.sample(Parental_Range_UMI, Offspring_SampleMutCount_RangeCount)
                    UMI_Select = UMI_Select + locals()['UMI_' + Ra]
                
                locals()['Parental_Mut_UMI' + Branch] = Parental_Mut_subset[Parental_Mut_subset['SampleName_UMI'].isin(UMI_Select)]
                Parental_UMI_MutCount = locals()['Parental_Mut_UMI' + Branch].value_counts('SampleName_UMI').reset_index(name = 'mut_count')[['mut_count']]
                Parental_UMI_MutCount['group'] = 'Parental'
                Parental_UMI_MutCount['Branch'] = 'Branch_' + Branch
                Parental_UMI_MutCount['num'] = i
                Parental_UMI_MutCount.to_csv(UMI_MutCount, sep = '\t', index=False, header=False)
            
            Parental_Mut_Branch_1 = locals()['Parental_Mut_UMI1']
            Parental_Mut_Branch_2 = locals()['Parental_Mut_UMI2']
            Parental_Mut_Branch_3 = locals()['Parental_Mut_UMI3']
            Mut_InDiffBranch_Count = len(set(Parental_Mut_Branch_1['mut_info'].unique()) & set(Parental_Mut_Branch_2['mut_info'].unique()) & set(Parental_Mut_Branch_3['mut_info'].unique()))
            
            MutRandom.write('Branch_1_2_3' + '\t' + 'Random' + '\t' + str(Mut_InDiffBranch_Count) + '\n')
        MutRandom.write('Branch_1_2_3' + '\t' + 'Observation' + '\t' + str(MutOffspring_Mut_BranchCount) + '\n')
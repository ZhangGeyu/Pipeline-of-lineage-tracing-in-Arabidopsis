Sanger_Seq_Raw_Data = 'Plant1_Progeny_Sanger_seq.fasta'
Output_Path_with_Consensus_Sanger_Seq = '/path/to/your/folder'

records = list(SeqIO.parse(Sanger_Seq_Raw_Data, "fasta"))

input_df = pd.DataFrame({"id": [record.id for record in records],\
    "sequence": [str(record.seq) for record in records]})
input_df['sample'] = input_df['id'].str.split('_').str.get(0)

Sample_list = input_df['sample'].unique()
print('# of total files: ' + str(len(Sample_list)))

def reverse_complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_complement_seq = ''.join(complement_dict[base] for base in reversed(sequence))
    return reverse_complement_seq

def find_common_sequence(seq1, seq2, min_length=50):
    common_sequences = []

    for i in range(len(seq1) - min_length + 1):
        for j in range(min_length, min(len(seq1) - i + 1, len(seq2) + 1)):
            if seq1[i:i+j] in seq2:
                common_sequences.append(seq1[i:i+j])
    return common_sequences

with open((Output_Path_with_Consensus_Sanger_Seq + 'R_F_Overlap.txt'),'w') as R_F_Overlap:
    R_F_Overlap.write('Sample\tConsensus_seq\tR_seq\tF_reverse_seq\tOverlap\n')
    
    with open((Output_Path_with_Consensus_Sanger_Seq + 'Sanger_ConsensusSequence.fa'),'w') as ConsensusSequence_fa:
        for sample in Sample_list:
            input_df_subset = input_df[input_df['sample'] == sample]
            R_seq = input_df_subset[input_df_subset['id'].str.contains('M13R')]['sequence'].unique()[0]
            F_seq = input_df_subset[input_df_subset['id'].str.contains('M13F')]['sequence'].unique()[0]
            F_seq_reverse = reverse_complement(F_seq)

            common_seqs = find_common_sequence(R_seq, F_seq_reverse, min_length = 50)
            longest_common_seq = max(common_seqs, key = len, default = 'None')
            print(f"Longest Overlap Sequence: {longest_common_seq}")

            if longest_common_seq != 'None':
                Consensus_seq = R_seq.split(longest_common_seq)[0] + longest_common_seq + F_seq_reverse.split(longest_common_seq)[-1]
                R_F_Overlap.write(sample + '\t' + Consensus_seq + '\t' + R_seq + '\t' + F_seq_reverse + '\t' + longest_common_seq + '\n')
                ConsensusSequence_fa.write('>' + sample + '\n' + Consensus_seq + '\n')

            elif longest_common_seq == 'None':
                Consensus_seq = R_seq
                R_F_Overlap.write(sample + '\t' + Consensus_seq + '\t' + R_seq + '\t' + F_seq_reverse + '\t' + longest_common_seq + '\n')

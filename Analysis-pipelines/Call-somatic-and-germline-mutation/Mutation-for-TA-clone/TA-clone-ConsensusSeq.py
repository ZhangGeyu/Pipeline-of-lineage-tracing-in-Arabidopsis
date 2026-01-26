from Bio import SeqIO
import pandas as pd

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

path = '../' # Path for the Raw data: Plant1_Progeny_Sanger_seq.fasta

Progeny_Sanger_Raw = list(SeqIO.parse((path + 'Plant1_Progeny_Sanger_seq.fasta'), "fasta"))
SampleName_list = list({record.id[:-5] for record in Progeny_Sanger_Raw})

with open((path + 'Plant1_Progeny_merge_seq.fasta'),'w') as R_F_Overlap_fasta:    
    for sample in SampleName_list:
        matched_records = [record for record in Progeny_Sanger_Raw if sample in record.id]
        R_seq = str([record for record in matched_records if 'M13R' in record.id][0].seq)
        F_seq = str([record for record in matched_records if 'M13F' in record.id][0].seq)
        F_seq_reverse = reverse_complement(F_seq)

        common_seqs = find_common_sequence(R_seq, F_seq_reverse, min_length = 50)
        longest_common_seq = max(common_seqs, key = len, default = 'None')
        #print(f"Longest Overlap Sequence: {longest_common_seq}")

        if longest_common_seq != 'None':
            Consensus_seq = R_seq.split(longest_common_seq)[0] + longest_common_seq + F_seq_reverse.split(longest_common_seq)[-1]
            R_F_Overlap_fasta.write('>' + sample + '\n' + Consensus_seq + '\n')

        elif longest_common_seq == 'None':
            Consensus_seq = R_seq
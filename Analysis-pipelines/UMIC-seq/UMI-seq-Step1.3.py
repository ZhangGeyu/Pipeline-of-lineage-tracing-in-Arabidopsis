import os
from Bio import SeqIO
import pandas as pd

Raw_fq_data = '../Plant1_Leaf_UMIC_seq.fq'
path = '../' # path for the /Demultiplex_fa folder

Demultiplex_files = [f for f in os.listdir(path + 'Demultiplex_fa/') if f.endswith('.fasta')]

BC_reads_id = []
BC_reads_sequences = []
BC_id = []

for Demultiplex_file in Demultiplex_files:
    SampleBC = Demultiplex_file.split('_')[-1].split('.')[0]
    BC_fasta_file = (path + 'Demultiplex_fa/' + Demultiplex_file)
    for record in SeqIO.parse(BC_fasta_file, "fasta"):
        BC_reads_id.append(record.id)
        BC_reads_sequences.append(str(record.seq))
        BC_id.append(SampleBC)

BC_reads_df = pd.DataFrame({'Reads_ID': BC_reads_id, 'Sequence': BC_reads_sequences, 'BC_ID': BC_id})

print('Finish: Sample Barcode reads to dataframe')
BC_reads_df['info'] = BC_reads_df['Reads_ID'] + '_' + BC_reads_df['BC_ID']
BC_reads_df = BC_reads_df.drop_duplicates('info').reset_index(drop = True)

# clean.fq to dataframe

data = []
with open(Raw_fq_data, 'r') as file:
    lines = file.readlines()
    buffer = []
    for line in lines:
        buffer.append(line.strip())
        if len(buffer) == 4:
            try:
                read_id = buffer[0]
                seq = buffer[1]
                plus = buffer[2]
                qual = buffer[3]
                data.append([read_id, seq, plus, qual])
            except Exception as e:
                print(f"Skipping invalid entry: {e}")
            buffer = []

clean_fq = pd.DataFrame(data, columns=['read_id', 'seq', '+', 'qual'])
clean_fq['read_id_split'] = clean_fq['read_id'].str.split(' ').str.get(0).str[1:]
clean_fq = clean_fq.drop_duplicates('read_id_split').reset_index(drop = True)

print('Finish: fastq to dataframe')

# split the clean.fq by sample barcode

for BC_ID, group in BC_reads_df.groupby("BC_ID"):
    Reads_ID_list = group["Reads_ID"].tolist()
    output_file = f"{BC_ID}.fq"
    print(BC_ID)
    clean_fq_subset = clean_fq[clean_fq['read_id_split'].isin(Reads_ID_list)].reset_index(drop = True)

    with open((path + 'BC_split_fq/' + output_file), "w") as output_fq:
        for i in range(len(clean_fq_subset['read_id_split'])):
            output_fq.write(clean_fq_subset['read_id'][i] + '\n')
            output_fq.write(clean_fq_subset['seq'][i] + '\n')
            output_fq.write(clean_fq_subset['+'][i] + '\n')
            output_fq.write(clean_fq_subset['qual'][i] + '\n')

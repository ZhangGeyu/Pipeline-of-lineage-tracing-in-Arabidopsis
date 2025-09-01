Sequence_Data_file = 'Plant1_Leaf_UMIC_seq.fq'
Barcode_Probe_file = 'barcode.probe.fasta'
Output_Reads_With_Different_BC = '/path/to/your/folder'
Output_Fastq_with_Different_BC = '/path/to/your/folder'

### Extraction of reads with sample barcodes by the barcode probe: 
# BC01 with --umi_len 15
'''
python UMIC-seq.py UMIextract --input Sequence_Data_file --probe Barcode_Probe_file --umi_loc up --umi_len 10 --output ExtractedBC.fasta
'''

### Split the Sequence_Data_file by sample barcode 

'''
python UMIC-seq_helper_guo.py demultiplex --barcodes barcodes.fasta --input ExtractedBC.fasta --output Output_Reads_With_Different_BC
'''

### Extract reads from Sequence_Data_file to generate .fastq file for each sample 

Demultiplex_files = [f for f in os.listdir(Output_Reads_With_Different_BC) if f.endswith('.fasta')]

BC_reads_id = []
BC_reads_sequences = []
BC_id = []

for Demultiplex_file in Demultiplex_files:
    SampleBC = Demultiplex_file.split('.')[0]
    BC_fasta_file = (Output_Reads_With_Different_BC + Demultiplex_file)
    for record in SeqIO.parse(BC_fasta_file, "fasta"):
        BC_reads_id.append(record.id)
        BC_reads_sequences.append(str(record.seq))
        BC_id.append(SampleBC)

BC_reads_df = pd.DataFrame({'Reads_ID': BC_reads_id, 'Sequence': BC_reads_sequences, 'BC_ID': BC_id})

print('Finish: Sample Barcode reads to dataframe')
BC_reads_df['info'] = BC_reads_df['Reads_ID'] + '_' + BC_reads_df['BC_ID']
BC_reads_df = BC_reads_df.drop_duplicates('info').reset_index(drop = True)

# Sequence_Data_file to dataframe

data = []
with open(Sequence_Data_file, 'r') as file:
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

# split the Sequence_Data_file by sample BC

for BC_ID, group in BC_reads_df.groupby("BC_ID"):
    Reads_ID_list = group["Reads_ID"].tolist()
    output_file = f"{BC_ID}.fq"
    print(BC_ID)
    clean_fq_subset = clean_fq[clean_fq['read_id_split'].isin(Reads_ID_list)].reset_index(drop = True)

    with open(Output_Fastq_with_Different_BC + output_file, "w") as output_fq:
        for i in range(len(clean_fq_subset['read_id_split'])):
            output_fq.write(clean_fq_subset['read_id'][i] + '\n')
            output_fq.write(clean_fq_subset['seq'][i] + '\n')
            output_fq.write(clean_fq_subset['+'][i] + '\n')
            output_fq.write(clean_fq_subset['qual'][i] + '\n')
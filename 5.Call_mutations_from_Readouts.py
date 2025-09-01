Output_File_for_Cluster_Full = '/path/to/your/folder'
Reads_Merge_File = '/path/to/your/folder/BC_UMI_merge.fasta'
reference_file = '/path/to/your/folder/reference.fa'

## Append the sample and readout numbers to the end of the read names and then merge them

input_files_name = [f.split('_')[1] for f in os.listdir(Output_File_for_Cluster_Full) if f.endswith('.pdf')]

reads_id = []
reads_sequences = []

for file_name in input_files_name:
    cluster_files_name = [n.split('.')[0].split('_')[1] for n in os.listdir(Output_File_for_Cluster_Full + 'UMIclusterfull_' + file_name) if n.endswith('.fasta')]
    print(file_name)
    for cluster in cluster_files_name:
        print(cluster)
        fasta_file = (Output_File_for_Cluster_Full + 'UMIclusterfull_' + file_name + '/cluster_' + cluster + '.fasta')
        for record in SeqIO.parse(fasta_file, "fasta"):
            reads_id.append((record.id + '-' + file_name + '-' + cluster))
            reads_sequences.append(str(record.seq))

reads_df = pd.DataFrame({'Reads_ID': reads_id, 'Sequence': reads_sequences})

with open(Reads_Merge_File, 'w') as fasta_file_merge:
    for index, row in reads_df.iterrows():
        reads_id = row['Reads_ID']
        seq = row['Sequence']
        fasta_file_merge.write(f'>{reads_id}\n{seq}\n')

## Mapping the reads to the reference.fa
'''
nohup minimap2 -ax map-ont -t 50 reference.fa BC_UMI_merge.fasta -o BC_UMI_merge.fasta.bam &
'''
## Sort the .bam file
'''
nohup samtools sort BC_UMI_merge.fasta.bam -o BC_UMI_merge.fasta.sorted.bam &
'''
## Split the bam file to single readout by the sample and readout numbers to the end of the read names

Sorted_Bam = '/path/to/your/folder/BC_UMI_merge.fasta.sorted.bam'
Output_Folder_for_Splited_Bam = '/path/to/your/folder'

bam_file = pysam.AlignmentFile(Sorted_Bam, 'rb')
header = str(bam_file.header)

SampleBC_Cluster_dict = {}

for read in bam_file:
    BC_cluster = '_'.join(str(read.query_name).split('-')[-2:])
    if BC_cluster in SampleBC_Cluster_dict:
        SampleBC_Cluster_dict[BC_cluster].append(read)
    else:
        SampleBC_Cluster_dict[BC_cluster] = [read]

SampleBC_Cluster_list = []
for SampleBC_Cluster, reads in SampleBC_Cluster_dict.items():
    SampleBC_Cluster_list.append(SampleBC_Cluster)
    print(SampleBC_Cluster)
    with open((Output_Folder_for_Splited_Bam + f"{SampleBC_Cluster}.bam"), "w") as file:
        file.write(header)
        for read in reads:
            file.write(f"{read.tostring()}\n")

## Command for batch processing to call mutation for all bam

reference_file = '/path/to/your/folder/reference.fa'
Output_Command_List_for_Mpileup = '/path/to/your/folder'

Sample_BC_Cluster_files = [f for f in os.listdir(Output_Folder_for_Splited_Bam) if f.endswith('.bam')]

print(len(Sample_BC_Cluster_files))
mpileup_command_list = []

for input_file in Sample_BC_Cluster_files:
    input_file_path = os.path.join(Output_Folder_for_Splited_Bam, input_file)
    output_file = input_file.replace(".bam", ".mpileup")
    output_file_path = os.path.join(Output_Folder_for_Splited_Bam, output_file)

    mpileup_command = f"samtools mpileup --max-depth 0 --output-BP --reference {reference_file} {input_file_path} -o {output_file_path}"
    mpileup_command_list.append(mpileup_command)

with open(Output_Command_List_for_Mpileup, "w") as file:
    for mpileup_command in mpileup_command_list:
        file.write("'" + mpileup_command + "'" + '\n')
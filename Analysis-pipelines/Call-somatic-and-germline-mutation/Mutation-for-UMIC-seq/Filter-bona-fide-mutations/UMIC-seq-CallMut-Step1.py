import os
from Bio import SeqIO

path = '../UMIclusterfull/'

input_files_name = [f.split('_')[1] for f in os.listdir(path) if f.endswith('.pdf')]

reads_id = []
reads_sequences = []

for file_name in input_files_name:
    cluster_files_name = [n.split('.')[0].split('_')[1] for n in os.listdir(path + 'UMIclusterfull_' + file_name) if n.endswith('.fasta')]
    print(file_name)
    for cluster in cluster_files_name:
        #print(cluster)
        fasta_file = (path + 'UMIclusterfull_' + file_name + '/cluster_' + cluster + '.fasta')
        for record in SeqIO.parse(fasta_file, "fasta"):
            reads_id.append((record.id + '-' + file_name + '-' + cluster))
            reads_sequences.append(str(record.seq))

reads_df = pd.DataFrame({'Reads_ID': reads_id, 'Sequence': reads_sequences})

with open((path + 'BC_UMI_merge.fasta'), 'w') as fasta_file_merge:
    for index, row in reads_df.iterrows():
        reads_id = row['Reads_ID']
        seq = row['Sequence']
        fasta_file_merge.write(f'>{reads_id}\n{seq}\n')
import pysam

path = '../' # Path to the folder containing BC_UMI_merge.fasta.sorted.bam

bam_file = pysam.AlignmentFile((path + 'BC_UMI_merge.fasta.sorted.bam'), 'rb')
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
    #print(SampleBC_Cluster)
    with open((path + 'BAM_and_mpileup_split/' + f"{SampleBC_Cluster}.bam"), "w") as file:
        file.write(header)
        for read in reads:
            file.write(f"{read.tostring()}\n")
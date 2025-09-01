Output_File_for_Extract_UMI = '/path/to/your/folder'
Output_Fastq_with_Different_BC = '/path/to/your/folder'
Output_File_for_Cluster_Full = '/path/to/your/folder'
Output_Command_List_for_Cluster_Full = '/path/to/your/folder'

### Command for batch processing to perform cluster full for each ExtractedUMIs.fasta: 

input_folder_1 = Output_File_for_Extract_UMI
input_folder_2 = Output_Fastq_with_Different_BC
output_folder = Output_File_for_Cluster_Full

input_files_name = [f.split('_')[1].split('.')[0] for f in os.listdir(input_folder_1) if f.endswith('.fasta')]
input_files_name.sort()
command_list = []

for input_file in input_files_name:
    input_file_path_1 = input_folder_1 + 'ExtractedUMIs_' + input_file + '.fasta'
    input_file_path_2 = input_folder_2 + input_file + '.fq'

    output_file = 'UMIclusterfull_' + input_file
    output_file_path = output_folder + output_file

    mpileup_command = f"'python UMIC-seq.py clusterfull --input {input_file_path_1} --reads {input_file_path_2} --aln_thresh 50 --size_thresh 10 --output {output_file_path} --stop_thresh 0'"
    command_list.append(mpileup_command)

with open(Output_Command_List_for_Cluster_Full, "w") as file:
    for mpileup_command in command_list:
        file.write(mpileup_command+'\n')
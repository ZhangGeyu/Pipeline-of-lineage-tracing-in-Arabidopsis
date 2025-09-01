Output_File_for_Extract_UMI = '/path/to/your/folder'
Output_File_for_Cluster_Test = '/path/to/your/folder'
Output_Command_List_for_Cluster_Test = '/path/to/your/folder'

###  Command for batch processing to perform clustering approximation for each ExtractedUMIs.fasta: 

input_folder = Output_File_for_Extract_UMI
output_folder = Output_File_for_Cluster_Test

input_files = [f for f in os.listdir(input_folder) if f.endswith('.fasta')]

command_list = []

for input_file in input_files:
    input_file_path = os.path.join(input_folder, input_file)
    output_file = 'UMIclustertest_' + input_file.split('.')[0].split('_')[1]
    output_file_path = os.path.join(output_folder, output_file)

    mpileup_command = f"'python UMIC-seq.py clustertest --input {input_file_path} --steps 20 60 5 --output {output_file_path}'"
    command_list.append(mpileup_command)

with open(Output_Command_List_for_Cluster_Test, "w") as file:
    for mpileup_command in command_list:
        file.write(mpileup_command+'\n')

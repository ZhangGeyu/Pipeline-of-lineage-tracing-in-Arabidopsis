Fastq_File_for_Each_Sample = '/path/to/your/folder'
Output_File_for_Extract_UMI = '/path/to/your/folder'
Output_Command_List_for_Extract_UMIs = '/path/to/your/folder'
probe_path = 'UMI_probe.fasta'

### Command for batch processing to identify UMI-containing reads in all samples using UMI probes

input_folder = Fastq_File_for_Each_Sample
output_folder = Output_File_for_Extract_UMI

input_files = [f for f in os.listdir(input_folder) if f.endswith('.fq')]

command_list = []

for input_file in input_files:
    input_file_path = os.path.join(input_folder, input_file)
    output_file = 'ExtractedUMIs_' + input_file.replace(".fq", ".fasta")
    output_file_path = os.path.join(output_folder, output_file)

    mpileup_command = f"'python UMIC-seq.py UMIextract --input {input_file_path} --probe {probe_path} --umi_loc down --umi_len 40 --output {output_file_path}'"
    command_list.append(mpileup_command)

with open(Output_Command_List_for_Extract_UMIs, "w") as file:
    for mpileup_command in command_list:
        file.write(mpileup_command+'\n')
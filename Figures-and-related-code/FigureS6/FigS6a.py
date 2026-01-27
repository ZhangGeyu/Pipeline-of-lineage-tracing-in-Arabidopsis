import re

def parse_vcf_file(input_file, output_file):

    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        f_out.write("Chr\tpos\tmutation\tmut_reads\ttotal_reads\tfrequency\n")
        for line in f_in:
            if line.startswith('##'):  # Skip header line 
                continue
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue
            #  Extract basic information
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]

            try:
                qual = float(fields[5])  # QUAL score
            except ValueError:
                continue
            if len(ref) != 1 or len(alt) != 1:  # Only process single nucleotide SNPs
                continue
            # Parse INFO field
            info_fields = fields[7].split(';')
            info_dict = {}
            for item in info_fields:
                if '=' in item:
                    key, value = item.split('=', 1)
                    info_dict[key] = value


            format_fields = fields[8].split(':')
            sample_data = fields[9].split(':')

  
            sample_dict = dict(zip(format_fields, sample_data))

            if 'AD' in sample_dict and 'DP' in sample_dict:
                try:
                    ad_values = sample_dict['AD'].split(',')
                    if len(ad_values) >= 2:
                        ref_reads = int(ad_values[0])
                        alt_reads = int(ad_values[1])
                        total_reads = int(sample_dict['DP'])
                        if alt_reads < 3:  # Filter: mutation reads must be at least 3
                            continue

                        qd = float(info_dict.get('QD', 0))
                        sor = float(info_dict.get('SOR', 999))  
                        fs = float(info_dict.get('FS', 999))   

                        # QUAL > 30, QD > 2, SOR < 3, FS < 60
                        if (qual > 30 and qd > 2 and sor < 3 and fs < 60):

                            mutation_freq = alt_reads / total_reads if total_reads > 0 else 0
                            mutation_type = f"{ref}>{alt}"
                            output_line = f"{chrom}\t{pos}\t{mutation_type}\t{alt_reads}\t{total_reads}\t{mutation_freq:.4f}\n"
                            f_out.write(output_line)
                except (ValueError, IndexError):
                    continue
def main():
    input_file = "UMIC_1.g.vcf"
    output_file = "UMIC_vcf_test_plant1.txt"

    parse_vcf_file(input_file, output_file)


if __name__ == "__main__":
    main()

               

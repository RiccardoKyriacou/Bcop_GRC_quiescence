import os
import argparse
from glob import glob
from collections import defaultdict

'''
This script
1) Takes TPM outputs from stringtie
2) Uses a dictionary to identify metadata tag with species/sex/tissue data
3) Outputs Scaffold TPM_value   Sample_name Germ/Soma
4) Looks at expression of GRC genes, comparing expression in soma, germ or both (gene_overlap_{development_stage}_TPM.tsv)
'''

# Usage : python3 get_TPM_values.py -t ../01_STAR_TPM_NM/outputs/

def process_stringtie_gtf(gtf_folder, output_file):
    # Find all stringtie GTF files in the specified folder
    gtf_files = glob(os.path.join(gtf_folder, "*TPM.gtf"))

    # Open the output file to write all entries
    with open(f"{output_file}.tsv", "w") as outf, open(f"{output_file}_GRC.tsv", "w") as outf2:
        print(f"Writing combined TPM to {output_file}")
        for gtf in gtf_files:
            # Extract sample ID (e.g., B1, B2, etc.) from the filename
            sample_id = os.path.basename(gtf).split('_')[0]

            # Parse the GTF file
            with open(gtf, 'r') as file:
                for line in file:
                    if line.startswith("#"):  # Skip headers
                        continue
                    fields = line.strip().split('\t')
                    scaffold_name = fields[0]
                    feature_type = fields[2]

                    # Skip unmapped scaffolds and exons
                    if "SCAFFOLD" in scaffold_name or "scaffold" in scaffold_name or feature_type == "exon":
                        continue

                    # Extract TPM and ref_gene_id values
                    values = fields[8]
                    TPM = None
                    ref_gene_id = None
                    attributes = values.strip().split(';')
                    for attr in attributes:
                        key_value = attr.strip().split(' ')
                        if len(key_value) < 2:  # check malformed
                            continue
                        key, value = key_value[0], key_value[1]
                        if key == "TPM":
                            TPM = float(value.strip('"'))
                        elif key == "ref_gene_id":
                            ref_gene_id = value.strip('"')
                        if TPM is not None and ref_gene_id is not None:
                            break
                    # Skip lines missing either value
                    if TPM is None or ref_gene_id is None:  
                        continue
                    # Write combined output file
                    outf.write(f"{scaffold_name}\t{TPM:.6f}\t{ref_gene_id}\t{sample_id}\n")

                    if "g" in ref_gene_id:
                        outf2.write(f"{scaffold_name}\t{TPM:.6f}\t{ref_gene_id}\t{sample_id}\n")

def main():
    parser = argparse.ArgumentParser(description="Process TPM values from stringtie outputs.")
    parser.add_argument("-t", "--TPM", type=str, help="Path to TPM GTF file folder", required=True)
    parser.add_argument("-o", "--output", type=str, default="combined_TPM_only.tsv", help="Path and name of output file")
    args = parser.parse_args()

    # parse out .tsv ending if added to -o 
    if ".tsv" in args.output:
        args.output = args.output.split(".")[0]

    # # Process the GTF files and write to the specified output file
    process_stringtie_gtf(args.TPM, args.output)

if __name__ == "__main__":
    main()
                                                                                                                                                              

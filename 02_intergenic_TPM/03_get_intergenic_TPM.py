import os
import argparse
from glob import glob
from collections import defaultdict

'''
This script
1) Takes TPM outputs from stringtie
2) Uses a dictionary to identify metadata tag with species/sex/tissue data
3) Outputs Scaffold TPM_value   Sample_name Germ/Soma
4) Looks at expression of GRC genes, comparing expression in soma, germ or both (gene_expression_overlap.tsv)
'''

# Usage : python3 03_get_intergenic_TPM.py -t outputs/

# Sample metadata dictionary
samples_info = {
    # Embryo 0-4h
    "ME1": {"species": "B_coprophila", "sex": "male", "tissue": "germ_embryo", "stage": "0-4h"},
    "ME2": {"species": "B_coprophila", "sex": "male", "tissue": "germ_embryo", "stage": "0-4h"},
    "ME3": {"species": "B_coprophila", "sex": "male", "tissue": "germ_embryo", "stage": "0-4h"},
    "FE1": {"species": "B_coprophila", "sex": "female", "tissue": "germ_embryo", "stage": "0-4h"},
    "FE2": {"species": "B_coprophila", "sex": "female", "tissue": "germ_embryo", "stage": "0-4h"},
    "FE3": {"species": "B_coprophila", "sex": "female", "tissue": "germ_embryo", "stage": "0-4h"},
    # Embryo 0-8h
    "ML1": {"species": "B_coprophila", "sex": "male", "tissue": "germ_embryo", "stage": "4-8h"},
    "ML2": {"species": "B_coprophila", "sex": "male", "tissue": "germ_embryo", "stage": "4-8h"},
    "ML3": {"species": "B_coprophila", "sex": "male", "tissue": "germ_embryo", "stage": "4-8h"},
    "FL1": {"species": "B_coprophila", "sex": "female", "tissue": "germ_embryo", "stage": "4-8h"},
    "FL2": {"species": "B_coprophila", "sex": "female", "tissue": "germ_embryo", "stage": "4-8h"},
    "FL3": {"species": "B_coprophila", "sex": "female", "tissue": "germ_embryo", "stage": "4-8h"},
    # Late larva/early pupa 
    "Fgerm1": {"species": "B_coprophila", "sex": "female", "tissue": "germ", "stage": "late-larva-early-pupa"},
    "Fgerm2": {"species": "B_coprophila", "sex": "female", "tissue": "germ", "stage": "late-larva-early-pupa"},
    "Fgerm3": {"species": "B_coprophila", "sex": "female", "tissue": "germ", "stage": "late-larva-early-pupa"},
    "Fbody1": {"species": "B_coprophila", "sex": "female", "tissue": "soma_carcass", "stage": "late-larva-early-pupa"},
    "Fbody2": {"species": "B_coprophila", "sex": "female", "tissue": "soma_carcass", "stage": "late-larva-early-pupa"},
    "Fbody3": {"species": "B_coprophila", "sex": "female", "tissue": "soma_carcass", "stage": "late-larva-early-pupa"},
    "Mgerm1": {"species": "B_coprophila", "sex": "male",   "tissue": "germ", "stage": "late-larva-early-pupa"},
    "Mgerm2": {"species": "B_coprophila", "sex": "male",   "tissue": "germ", "stage": "late-larva-early-pupa"},
    "Mgerm3": {"species": "B_coprophila", "sex": "male",   "tissue": "germ", "stage": "late-larva-early-pupa"},
    "Mbody1": {"species": "B_coprophila", "sex": "male",   "tissue": "soma_carcass", "stage": "late-larva-early-pupa"},
    "Mbody2": {"species": "B_coprophila", "sex": "male",   "tissue": "soma_carcass", "stage": "late-larva-early-pupa"},
    "Mbody3": {"species": "B_coprophila", "sex": "male",   "tissue": "soma_carcass", "stage": "late-larva-early-pupa"},
    # Adult samples for B_coprophila
    "B1":  {"species": "B_coprophila", "sex": "male", "tissue": "germ_testes",  "stage": "adult"},
    "B2":  {"species": "B_coprophila", "sex": "male", "tissue": "germ_testes",  "stage": "adult"},
    "B3":  {"species": "B_coprophila", "sex": "male", "tissue": "germ_testes",  "stage": "adult"},
    "B7":  {"species": "B_coprophila", "sex": "male", "tissue": "soma_carcass", "stage": "adult"},
    "B8":  {"species": "B_coprophila", "sex": "male", "tissue": "soma_carcass", "stage": "adult"},
    "B9":  {"species": "B_coprophila", "sex": "male", "tissue": "soma_carcass", "stage": "adult"},
    "B10": {"species": "B_coprophila", "sex": "gynogenic_female", "tissue": "germ_ovaries", "stage": "adult"},
    "B11": {"species": "B_coprophila", "sex": "gynogenic_female", "tissue": "germ_ovaries", "stage": "adult"},
    "B12": {"species": "B_coprophila", "sex": "gynogenic_female", "tissue": "germ_ovaries", "stage": "adult"},
    "B16": {"species": "B_coprophila", "sex": "gynogenic_female", "tissue": "soma_carcass", "stage": "adult"},
    "B17": {"species": "B_coprophila", "sex": "gynogenic_female", "tissue": "soma_carcass", "stage": "adult"},
    "B18": {"species": "B_coprophila", "sex": "gynogenic_female", "tissue": "soma_carcass", "stage": "adult"},
    "B19": {"species": "B_coprophila", "sex": "androgenic_female", "tissue": "germ_ovaries", "stage": "adult"},
    "B20": {"species": "B_coprophila", "sex": "androgenic_female", "tissue": "germ_ovaries", "stage": "adult"},
    "B21": {"species": "B_coprophila", "sex": "androgenic_female", "tissue": "germ_ovaries", "stage": "adult"},
    "B25": {"species": "B_coprophila", "sex": "androgenic_female", "tissue": "soma_carcass", "stage": "adult"},
    "B26": {"species": "B_coprophila", "sex": "androgenic_female", "tissue": "soma_carcass", "stage": "adult"},
    "B27": {"species": "B_coprophila", "sex": "androgenic_female", "tissue": "soma_carcass", "stage": "adult"},
}


def process_stringtie_gtf(gtf_folder, samples_info, output_file):
    # Find all stringtie GTF files in the specified folder
    gtf_files = glob(os.path.join(gtf_folder, "*intergenic.gtf"))

    # Open the output file to write all entries
    with open(output_file, "w") as outf:
        print(f"Writing combined TPM to {output_file}")
        for gtf in gtf_files:
            # Extract sample ID (e.g., B1, B2, etc.) from the filename
            sample_id = os.path.basename(gtf).split('_')[0]

            # Get sample metadata
            sample_metadata = samples_info.get(sample_id)
            if not sample_metadata:
                continue  # Skip if metadata is missing
            species = sample_metadata['species']
            sex = sample_metadata['sex']
            tissue_type = sample_metadata['tissue']
            stage = sample_metadata['stage']

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
                    outf.write(f"{scaffold_name}\t{TPM:.6f}\t{ref_gene_id}\t{species}\t{sex}\t{tissue_type}\t{stage}\t{sample_id}\n")

def main():
    parser = argparse.ArgumentParser(description="Process TPM values from stringtie outputs.")
    parser.add_argument("-t", "--TPM", type=str, help="Path to TPM GTF file folder", required=True)
    parser.add_argument("-o", "--output", type=str, default="combined_intergenic_TPM.tsv", help="Path and name of output file")
    args = parser.parse_args()
    
    # Process the GTF files and write to the specified output file
    process_stringtie_gtf(args.TPM, samples_info, args.output)

if __name__ == "__main__":
    main()



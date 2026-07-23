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

# Sample metadata dictionary
samples_info = {
    # Embryo 0-4h
    "ME1": {"species": "B_coprophila", "sex": "male", "tissue": "germ_embryo", "stage": "0-4h"},
    "ME2": {"species": "B_coprophila", "sex": "male", "tissue": "germ_embryo", "stage": "0-4h"},
    "ME3": {"species": "B_coprophila", "sex": "male", "tissue": "germ_embryo", "stage": "0-4h"},
    "FE1": {"species": "B_coprophila", "sex": "female", "tissue": "germ_embryo", "stage": "0-4h"},
    "FE2": {"species": "B_coprophila", "sex": "female", "tissue": "germ_embryo", "stage": "0-4h"},
    "FE3": {"species": "B_coprophila", "sex": "female", "tissue": "germ_embryo", "stage": "0-4h"},
    # Embryo 4-8h
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
    gtf_files = glob(os.path.join(gtf_folder, "*TPM.gtf"))

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


def gene_overlap(combined_TPM_file, tpm_threshold, development_stage):
    """
    Process the combined TPM file, categorize genes based on expression,
    and return both a summary and the overlap data as a list of tuples.
    Each tuple contains (gene_id, tissue, sex).
    """
    # Dictionary to store gene expression data by tissue type and sex
    gene_expression = defaultdict(lambda: {"male": {"soma": 0, "germ": 0},
                                           "female": {"soma": 0, "germ": 0}})
    
    with open(combined_TPM_file) as f:
        for line in f:
            scaffold_name, TPM, ref_gene_id, species, sex, tissue_type, stage, sample_id = line.strip().split("\t")
            tpm = float(TPM)

            if stage == development_stage:
                # Extract GRC gene expression
                if "g" in ref_gene_id and tpm > tpm_threshold:
                    # Get sex key 
                    sex_key = "male" if sex == "male" else "female"
                    # Get tissue key 
                    tissue_key = "soma" if "soma" in tissue_type else "germ"
                    gene_expression[ref_gene_id][sex_key][tissue_key] += tpm

    # Categorise genes based on expression in tissue and sex
    categories = {
        "male_soma_only": [],
        "male_germ_only": [],
        "male_both": [],
        "female_soma_only": [],
        "female_germ_only": [],
        "female_both": [],
        "both-sexes_soma": [],
        "both-sexes_germ": [],
        "both-sexes_both": []
    }

    for gene_id, expression in gene_expression.items():
        male_soma = expression["male"]["soma"]
        male_germ = expression["male"]["germ"]
        female_soma = expression["female"]["soma"]
        female_germ = expression["female"]["germ"]

        # Male-specific categories
        if male_soma > 0 and male_germ == 0 and female_soma == 0 and female_germ == 0:
            categories["male_soma_only"].append(gene_id)
        elif male_germ > 0 and male_soma == 0 and female_soma == 0 and female_germ == 0:
            categories["male_germ_only"].append(gene_id)
        elif male_soma > 0 and male_germ > 0 and female_soma == 0 and female_germ == 0:
            categories["male_both"].append(gene_id)

        # Female-specific categories
        elif female_soma > 0 and female_germ == 0 and male_soma == 0 and male_germ == 0:
            categories["female_soma_only"].append(gene_id)
        elif female_germ > 0 and female_soma == 0 and male_soma == 0 and male_germ == 0:
            categories["female_germ_only"].append(gene_id)
        elif female_soma > 0 and female_germ > 0 and male_soma == 0 and male_germ == 0:
            categories["female_both"].append(gene_id)

        # Both sexes
        elif (male_soma > 0 or female_soma > 0) and (male_germ == 0 and female_germ == 0):
            categories["both-sexes_soma"].append(gene_id)
        elif (male_germ > 0 or female_germ > 0) and (male_soma == 0 and female_soma == 0):
            categories["both-sexes_germ"].append(gene_id)
        elif (male_soma > 0 or female_soma > 0) and (male_germ > 0 or female_germ > 0):
            categories["both-sexes_both"].append(gene_id)

    # Build the overlap data in memory (list of tuples)
    overlap_data = []
    for category, genes in categories.items():
        for gene_id in genes:
            sex = category.split("_")[0]
            tissue = category.split("_")[1]
            overlap_data.append((gene_id, tissue, sex))
    
    # Return both a summary and the overlap data (without writing an intermediate file)
    return overlap_data

def geneID_to_scaffold(combined_TPM_file):
    geneID_scaffold_dict = {}
    with open(combined_TPM_file) as f:
        for line in f:
            columns = line.strip().split("\t")
            scaffold = columns[0]
            gene_id = columns[2]
            geneID_scaffold_dict[gene_id] = scaffold
    
    return geneID_scaffold_dict

def add_TPM(overlap_data, geneID_scaffold_dict, combined_f, output_f, development_stage):
    """
    Combine the overlap data with TPM values from the combined file and write the final output.
    """
    # Dictionary to store gene_id with sample_id:TPM mappings
    gene_TPM_dict = defaultdict(list)
    
    # Read the combined file and store TPM values by gene_id
    with open(combined_f) as f1:
        for line in f1:
            scaffold_name, TPM, ref_gene_id, species, sex, tissue_type, stage, sample_id = line.strip().split("\t")
            # Only keep for same development stage
            if stage == development_stage:
                gene_TPM_dict[ref_gene_id].append(f"{sample_id}:{TPM}")
    
    # Write the final output file
    with open(output_f, "w") as outf:
        for gene_id, tissue, sex in overlap_data:
            sample_TPM_str = "\t".join(gene_TPM_dict.get(gene_id, []))
            scaffold = geneID_scaffold_dict.get(gene_id, [])
            outf.write(f"{gene_id}\t{scaffold}\t{tissue}\t{sex}\t{sample_TPM_str}\n")

def main():
    parser = argparse.ArgumentParser(description="Process TPM values from stringtie outputs.")
    parser.add_argument("-t", "--TPM", type=str, help="Path to TPM GTF file folder", required=True)
    parser.add_argument("-o", "--output", type=str, default="combined_TPM_only.tsv", help="Path and name of output file")
    args = parser.parse_args()
    
    # # Process the GTF files and write to the specified output file
    process_stringtie_gtf(args.TPM, samples_info, args.output)

    geneID_scaffold_dict = geneID_to_scaffold(args.output)
    
    # List of development development_stages to process
    development_stages = ["0-4h", "4-8h", "late-larva-early-pupa", "adult"]
    for stage in development_stages:
        # Get overlap data in memory (no intermediate file is written)
        overlap_data = gene_overlap(args.output, tpm_threshold=0, development_stage=stage)
        # Write the final output file combining overlap and TPM values
        final_filename = f"gene_overlap_{stage}_TPM.tsv"
        add_TPM(overlap_data, geneID_scaffold_dict, args.output, final_filename, development_stage=stage)

if __name__ == "__main__":
    main()
                                                                                                                                                              
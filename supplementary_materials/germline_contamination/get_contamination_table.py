from Bio import SeqIO
from glob import glob
import pandas as pd
import argparse
import os

"""
"""

# Get gene IDs for gernm-specifc non-GRC genes
def get_gene_ids(gene_file):
    ID_name_dict = {}
    with open(gene_file) as f:
        for line in f:
            homologue_name, gene_ID = line.strip().split("\t")
            ID_name_dict[gene_ID] = homologue_name
    return ID_name_dict

# Look at combined_TPM_only.tsv (generated in step 01_STAR_TPM)
def check_expression(ID_name_dict, TPM_file):
    with open(TPM_file) as f, open("contamination_table.tsv", "w") as outf:
        for line in f:
            scaffold, TPM, ID, sp, sex, tissue, stage, sample_ID = line.strip().split("\t")

            if float(TPM) > float(0.222534):
                if ("germ" in tissue) and ("female" in sex):
                    tissue = "germ_ovaries"
                    sex = "female"
                elif ("germ" in tissue) and ("male" in sex):
                    tissue = "germ_testes"

                if ID in ID_name_dict:
                    homologue_name = ID_name_dict[ID]
                    outf.write(f"{homologue_name}\t{ID}\t{tissue}\t{TPM}\t{sex}\t{stage}\n")
            else:
                continue

def main():
    parser = argparse.ArgumentParser(description="Get expressed genes for BLAST")
    parser.add_argument("-g", "--gene", type=str, required=True, help="Path to list of genes")
    parser.add_argument("-t", "--tpm", type=str, required=True, help="Path to combined_TPM_only.tsv (generated in step 01_STAR_TPM)")
    args = parser.parse_args()

    ID_name_dict = get_gene_ids(args.gene)
    check_expression(ID_name_dict, args.tpm)

if __name__ == "__main__":
    main()

from Bio import SeqIO
from glob import glob
import pandas as pd
import argparse
import os

"""
This script extracts gene sequences from a genome FASTA file using annotations provided in a GFF3 file.
It is designed to retrieve either all genes or a specified subset and outputs their full genomic sequences
in FASTA format. This is particularly useful for preparing sequences for downstream analysis such as BLAST.

Features:
- Parses a GFF3 file to identify gene features.
- Extracts corresponding DNA sequences from a genome FASTA file.
- Supports both forward and reverse strands.
- Accepts an optional gene list to extract only specific genes.
- Outputs a FASTA file containing sequences of each gene.

Usage:
    python 01_get_expressed_genes_for_BLAST.py \
        -g path/to/genome.fa \
        -a path/to/annotations.gff3 \
        [-l path/to/gene_list.txt] \
        [-o output_prefix]

Arguments:
    -g, --genome       (Required) Path to the genome FASTA file
    -a, --annotation   (Required) Path to the GFF3 annotation file
    -l, --list         (Optional) File with gene IDs to extract (one per line)
    -o, --output       (Optional) Output FASTA file prefix (default: annotated_genes)

Output:
    A FASTA file named <output_prefix>.fasta containing the extracted gene sequences
"""

# 1) Optional: parse GFF file to extract gene IDs
def get_GFF_gene_IDs(gff_file):
    gene_ids = set()
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            feature_type = fields[2]
            attributes = fields[8]
            if feature_type == 'gene':
                for attribute in attributes.split(';'):
                    if attribute.startswith('ID='):
                        gene_id = attribute.split('=')[1]
                        gene_ids.add(gene_id)
    return gene_ids 

# Function to extract and stitch exons for specified GRC genes
def parse_gff(genome_fasta, gff_file, gene_lst):
    GRC_gene_seqs = {}
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

    with open(gff_file, 'r') as gff:
        for line in gff:
            if line.startswith("#"):
                continue

            columns = line.strip().split('\t')
            if len(columns) < 9:
                continue

            scaffold = columns[0]
            feature_type = columns[2]
            start = int(columns[3])
            end = int(columns[4])
            strand = columns[6]
            attributes = columns[8]

            if feature_type == "gene":
                gene_ID = attributes.split("=")[1].split(";")[0]
                if gene_ID in gene_lst:
                    gene_sequence = genome[scaffold][start-1:end]
                    if strand == "-":
                        gene_sequence = gene_sequence.reverse_complement()
                    GRC_gene_seqs[gene_ID] = gene_sequence
    return GRC_gene_seqs

# Function to write to FASTA
def combined_fasta_output(output_file, gene_sequences):
    with open(f"{output_file}.fasta", 'w') as out_fasta:
        for gene_id, seq in gene_sequences.items():
            out_fasta.write(f">{gene_id}\n{str(seq.seq)}\n")

# Main block
def main():
    parser = argparse.ArgumentParser(description="Get expressed genes for BLAST")
    parser.add_argument("-g", "--genome", type=str, required=True, help="Path to genome FASTA")
    parser.add_argument("-a", "--annotation", type=str, required=True, help="Path to GFF3 file")
    parser.add_argument("-l", "--list", type=str, help="Optional: file with list of gene IDs")
    parser.add_argument("-o", "--output", type=str, default="annotated_genes", help="Output FASTA file prefix (default: annotated_genes)")
    args = parser.parse_args()

    if args.list:
        with open(args.list) as f:
            gene_lst = [line.strip() for line in f if line.strip()]
        custom_gene_sequences = parse_gff(args.genome, args.annotation, gene_lst)
        combined_fasta_output(args.output, custom_gene_sequences)
    else:
        gene_lst = get_GFF_gene_IDs(args.annotation)
        all_gene_sequences = parse_gff(args.genome, args.annotation, gene_lst)
        combined_fasta_output(args.output, all_gene_sequences)


if __name__ == "__main__":
    main()

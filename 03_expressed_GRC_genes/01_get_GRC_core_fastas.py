from Bio import SeqIO
from glob import glob
import pandas as pd
import argparse
import os

'''
This script
1) Using gene ID's parses out the gene sequences from B_cop GRC genome
2) Outputs three FASTA files for GRC genes expressed in soma, germ, and both
3) Parses out all core genes from annotated B_cop genome
4) Produces combined fasta file contining all core chromosome genes
'''

# 1) Get list of GRC genes
def get_expressed_gene_ids(expressed_GRC_genes_summary):
    GRC_gene_lst = [] 
    with open(expressed_GRC_genes_summary, "r") as f:
        next(f)  # skip header line
        for line in f:
            columns = line.strip().split("\t")
            gene_id = columns[0]
            GRC_gene_lst.append(gene_id)
    return set(GRC_gene_lst)       
         
# 2) Get BLAST quieries for GRC genes and all core genes   
# Function to extract and stitch exons for specified GRC genes
def parse_gff(genome_fasta, gff_file, expressed_GRC_gene_lst):
    # Dictionary to store stitched exon sequences by gene
    GRC_gene_seqs = {}

    # Load the genome sequences into a dictionary
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

    # Parse GFF file and extract gene sequences
    with open(gff_file, 'r') as gff:
        for line in gff:
            # Skip any headers/comments
            if line.startswith("#"):
                continue
            
            columns = line.strip().split('\t')
            scaffold = columns[0] 
            feature_type = columns[2]  # Feature type (gene, exon, etc.)
            start = int(columns[3])  # Start position (1-based indexing)
            end = int(columns[4])  # End position
            strand = columns[6]  # Strand (+ or -)

            # Only consider gene
            if feature_type == "gene":
                # Extract the parent gene ID
                gene_ID = columns[8].split("=")[1].split(";")[0]  # Parse ID

                if gene_ID in expressed_GRC_gene_lst:
                    # Get the sequence from the genome
                    gene_sequence = genome[scaffold][start-1:end]  # Adjust to 0-based indexing

                    # Reverse complement if on the negative strand
                    if strand == "-":
                        gene_sequence = gene_sequence.reverse_complement()

                    # Append exonic sequence to the corresponding gene
                    GRC_gene_seqs[gene_ID] = gene_sequence

    return GRC_gene_seqs

# Function to write transcripts to FASTA format
def combined_fasta_output(output_file, transcript_sequences):
    with open(f"{output_file}.fasta", 'w') as out_fasta:
        for gene_id, seq in transcript_sequences.items():
            # Write each transcript in FASTA format
            out_fasta.write(f">{gene_id}\n{str(seq.seq)}\n")

# 2) Functions to generate FASTA file of all core-genome genes
def parse_core_genes_from_gff(genome_fasta, gff_file):
    # Dictionary to store stitched exon sequences by gene
    geneID_sequences = {}

    # Load the genome sequences into a dictionary
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

    # Parse GFF file and extract exonic sequences
    with open(gff_file, 'r') as gff:
        for line in gff:
            # Skip any headers/comments
            if line.startswith("#"):
                continue
            
            columns = line.strip().split('\t')
            scaffold = columns[0] 
            feature_type = columns[2]  # Feature type (gene, exon, etc.)
            start = int(columns[3])  # Start position (1-based indexing)
            end = int(columns[4])  # End position
            strand = columns[6]  # Strand (+ or -)

            # Only consider full genes
            if feature_type == "gene":
                # Extract the parent gene ID
                gene_ID = columns[8].split("=")[1].split(";")[0]  # Parse ID

                if "s" in gene_ID:
                    # Get the sequence from the genome
                    gene_sequence = genome[scaffold][start-1:end]  # Adjust to 0-based indexing

                    # Reverse complement if on the negative strand
                    if strand == "-":
                        gene_sequence = gene_sequence.reverse_complement()

                    geneID_sequences[gene_ID] = gene_sequence
                        
    return geneID_sequences

def main():
    parser = argparse.ArgumentParser(description="get expressed genes for BLAST")
    parser.add_argument("-g", "--genome", type=str, help="Path to genome to parse", required=True)
    parser.add_argument("-a", "--annotation", type=str, help="Path to annotation", required=True)
    parser.add_argument("-s", "--summary", type=str, help="Path to Expressed_GRC_genes_summary.tsv", required=True)
    args = parser.parse_args()

    # 1) Expression table and GRC gene sequences
    expressed_GRC_genes = get_expressed_gene_ids(args.summary)
    GRC_gene_seqs = parse_gff(args.genome, args.annotation, expressed_GRC_genes)
    combined_fasta_output("GRC_genes_to_BLAST", GRC_gene_seqs)

    # 2) Core GRC gene fasta
    core_gene_sequences = parse_core_genes_from_gff(args.genome, args.annotation)
    combined_fasta_output("core_genes", core_gene_sequences)


if __name__ == "__main__":
    main()

# TODO add BLAST step 
# python 01_get_GRC_core_fastas.py -g ../../../Annotations/idLycInge5.1.primary.masked.fa -a ../../../Annotations/ling_core_GRC.gff -s ../outputs/Expressed_GRC_genes_summary.tsv

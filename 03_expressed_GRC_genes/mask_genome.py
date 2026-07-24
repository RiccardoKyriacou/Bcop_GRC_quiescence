import argparse
from Bio import SeqIO

"""
Script to mask scaffold from a genome.
Default: Masks GRC genome.
Toggle: use -m core to mask out core genome instead.
"""

def mask_genome_region(genome_fasta, outfile, region_to_mask):
    # Load the genome sequences into a dictionary
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

    # Dictionary comprehension to store scaffold sequences as lists of nucleotides
    genome_seqs = {scaffold.id: list(scaffold.seq) for scaffold in genome.values()}

    for scaffold_id, record in genome.items():
        if region_to_mask == "GRC" and "GRC" in scaffold_id:  # Mask GRC scaffolds
            genome_seqs[scaffold_id] = ['N'] * len(record.seq)
        elif region_to_mask == "core" and "GRC" not in scaffold_id:  # Mask core scaffolds
            genome_seqs[scaffold_id] = ['N'] * len(record.seq)

    # Write out masked genome
    with open(outfile, "w") as f:
        for scaffold_id, seq_list in genome_seqs.items():
            masked_genome_seq = "".join(seq_list)
            f.write(f">{scaffold_id}\n{masked_genome_seq}\n")

def main():
    parser = argparse.ArgumentParser(description="Mask scaffolds in a genome based on the specified region.")
    parser.add_argument("-g", "--genome", type=str, help="Path to the genome file (FASTA format)", required=True)
    parser.add_argument("-o", "--output", type=str, help="Output file name", required=True)
    parser.add_argument("-m", "--mask", type=str, choices=["core", "GRC"], help="Region to mask: core or GRC", required=True)
    args = parser.parse_args()

    # Ensure the output filename ends with .fasta
    if not args.output.endswith(".fasta"):
        args.output += ".fasta"

    # Mask the genome based on the selected region
    mask_genome_region(args.genome, args.output, args.mask)

if __name__ == "__main__":
    main()

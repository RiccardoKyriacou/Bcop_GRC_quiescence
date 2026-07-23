import argparse
from Bio import SeqIO
from Bio.Seq import Seq

'''
Script to generate translated transcripts (proteins) when given a list of genes, genome, and annotation.
This version applies phase correction after splicing rather than on each CDS segment.
'''

# Function to get gene IDs for BLAST candidates remains the same
def get_expressed_gene_ids(GRC_expression_file):
    GRC_gene_lst = [] 
    with open(GRC_expression_file) as f:
        for line in f:
            if line.startswith("gene_id"):
                pass
            else:
                columns = line.strip().split("\t")
                gene_id = columns[0]
                GRC_gene_lst.append(gene_id)
    return set(GRC_gene_lst)

# Updated function to parse the GFF file and generate spliced CDS sequences
def parse_gff(genome_fasta, gff_file, GRC_gene_lst):
    """
    This function parses a GFF file, extracts only the CDS features, and groups them by transcript.
    It then splices the segments together and applies phase correction as a final trimming step.
    
    For positive-strand transcripts:
      - CDS features are sorted in ascending order.
      - The frame of the first CDS is used to trim the beginning of the spliced sequence.
    
    For negative-strand transcripts:
      - CDS features are sorted in ascending order.
      - The biologically first CDS is the one with the highest coordinate (last in sorted order).
      - After concatenation, the entire sequence is reverse complemented and then trimmed by the phase offset of that last (biologically first) CDS.
    """
    # Load the genome sequences into a dictionary
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    
    # Dictionary to store CDS features grouped by transcript
    cds_dict = {}
    
    with open(gff_file, 'r') as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            
            columns = line.strip().split('\t')
            scaffold = columns[0]
            feature_type = columns[2]
            start = int(columns[3])     # 1-based start
            end = int(columns[4])       # End coordinate
            strand = columns[6]         # Strand information: + or -
            phase = columns[7]
            phase = int(phase) if phase != "." else 0
            attributes = columns[8]
            
            # Only process CDS features
            if feature_type == "CDS":
                # Extract transcript ID from Parent attribute
                # Example: "ID=g2596.t1.CDS1;Parent=g2596.t1;"
                transcript_id = attributes.split("Parent=")[1].split(";")[0]
                parent_gene_id = transcript_id.split(".")[0]
                if parent_gene_id in GRC_gene_lst:
                    if transcript_id not in cds_dict:
                        cds_dict[transcript_id] = {"strand": strand, "scaffold": scaffold, "features": []}
                    cds_dict[transcript_id]["features"].append((start, end, phase))
    
    transcript_sequences = {}
    for transcript_id, info in cds_dict.items():
        strand = info["strand"]
        scaffold = info["scaffold"]
        features = info["features"]
        
        # For phase correction treat positive and negative strands differently.
        # Sort features in ascending order based on their genomic coordinate.
        features.sort(key=lambda x: x[0]) 
        
        # Build full spliced sequence by concatenating the full nucleotide segments.
        spliced_seq = ""
        for start, end, phase in features:
            # Extract full segment (using 0-based indexing)
            seg = genome[scaffold].seq[start - 1 : end]
            spliced_seq += str(seg)
        
        # Convert concatenated sequence to a Biopython Seq object.
        spliced_seq = Seq(spliced_seq)
        
        if strand == "+":
            # For positive strand, the biologically first CDS is the first in sorted order.
            offset = features[0][2]  # phase from the first CDS feature
            # Trim the beginning of the spliced sequence by the offset.
            spliced_seq = spliced_seq[offset:]
        else:
            # For negative strand, the biologically first CDS is the last one in genomic order.
            offset = features[-1][2]  # phase from the last CDS feature (which is biologically first)
            # Reverse complement the concatenated sequence.
            spliced_seq = spliced_seq.reverse_complement()
            # Trim the beginning of the reversed sequence by the offset.
            spliced_seq = spliced_seq[offset:]
        
        transcript_sequences[transcript_id] = spliced_seq
        
    return transcript_sequences

def translate_transcripts(transcript_sequences):
    """
    Translates the spliced nucleotide sequences into amino acid sequences.
    """
    aa_sequences = {}
    for transcript_id, nucl_seq in transcript_sequences.items():
        try:
            # Translate assuming standard codon table and that input is in-frame
            aa_seq = nucl_seq.translate(to_stop=False)  # Don't stop at first stop codon
        except Exception as e:
            print(f"Warning: {transcript_id} encountered an issue with translation: {e}.")
            aa_seq = "TRANSLATION_ERROR"
        aa_sequences[transcript_id] = aa_seq
    return aa_sequences

# Function to write the amino acid sequences to FASTA format
def combined_fasta_output(output_file, aa_sequences):
    with open(f"{output_file}.fasta", 'w') as out_fasta:
        for transcript_id, aa_seq in aa_sequences.items():
            out_fasta.write(f">{transcript_id}\n{str(aa_seq)}\n")

def main():
    parser = argparse.ArgumentParser(description="Extract, stitch and translate CDS sequences for GRC genes.")
    parser.add_argument("-t", "--table", type=str, help="Path to expression table", required=True)
    parser.add_argument("-g", "--genome", type=str, help="Path to genome fasta", required=True)
    parser.add_argument("-a", "--annotation", type=str, help="Path to GFF annotation", required=True)
    args = parser.parse_args()

    GRC_gene_lst = get_expressed_gene_ids(args.table)
    transcript_sequences = parse_gff(args.genome, args.annotation, GRC_gene_lst)
    aa_sequences = translate_transcripts(transcript_sequences)
    combined_fasta_output("GRC_transcripts", aa_sequences)
    combined_fasta_output("GRC_genes", transcript_sequences)

if __name__ == "__main__":
    main()

# Usage:
# python 03_get_expressed_GRC_transcripts.py -t ../outputs/Corrected_expressed_GRC_genes.tsv -g ../../../../../Annotations/idBraCopr2.1.primary.masked.fa -a ../../../../../Annotations/bcop_core_GRC.gff3

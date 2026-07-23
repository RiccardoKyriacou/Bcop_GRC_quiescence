#!/usr/bin/env python3
import argparse
import pandas as pd

def parse_blast_file(filename):
    # Get best hit for each gene
    best_hit = {}
    with open(filename) as f:
        for line in f:
            columns = line.rstrip("\n").split("\t")
            qseqid = columns[0]
            sseqid = columns[1]
            try:
                pident = float(columns[2])
                length = int(columns[3])
                mismatch = int(columns[4])
                evalue = float(columns[10])
                qlen = int(columns[12])
            except ValueError:
                continue  # skip lines that don’t parse correctly
           
            # Calculate coverage (in percentage)
            coverage = round((length / qlen) * 100, 3)
            
            # Save this hit if either no previous hit exists or if it is better.
            # "Better" is defined as having a higher coverage, or if coverage is equal, a lower e-value.
            if qseqid not in best_hit:
                best_hit[qseqid] = {
                    'sseqid': sseqid,
                    'coverage': coverage,
                    'pident': pident,
                    'evalue': evalue,
                    'qlen': qlen,
                    'mismatch': mismatch}
            else:
                prev = best_hit[qseqid]
                if (coverage > prev['coverage'] or 
                    (coverage == prev['coverage'] and evalue < prev['evalue'])):
                    best_hit[qseqid] = {
                        'sseqid': sseqid,
                        'coverage': coverage,
                        'pident': pident,
                        'evalue': evalue,
                        'qlen': qlen,
                        'mismatch': mismatch}
                
    return best_hit

def get_BLAST_table(expression_table, gene_blast, genome_blast, outputf):
    # Read expression table and store info
    GRC_gene_info = []
    with open(expression_table) as f:
        next(f)  # Skip the first line (header)
        for line in f:
            line = line.rstrip("\n")
            columns = line.split("\t")
            GRC_gene_info.append(columns)

    # Parse both BLAST files into dictionaries keyed by gene ID.
    core_hits = parse_blast_file(gene_blast)
    genome_hits = parse_blast_file(genome_blast)

    # Write output table
    with open(outputf, "w") as outf:
        #Write table header
        outf.write(f"gene_id\tscaffold\ttissue_category\tdevelopment_stage\tsex\tgerm/soma\tmean_germ_TPM/mean_soma_TPM\tHit_id\tCoverage\t%Identity\tE-Value\tHit_length\tMismatches\n")
        for row in GRC_gene_info:

            gene_id = row[0]
            scaffold= row[1]
            tissue_expr= row[2]
            development_stage= row[3]
            sex= row[4]
            count_expr_germ= row[5]
            count_expr_soma= row[6]
            mean_germ_TPM= float(row[7])
            mean_germ_TPM = round(mean_germ_TPM, 2)
            mean_soma_TPM= float(row[8])
            mean_soma_TPM = round(mean_soma_TPM, 2)

            # Look first in the core BLAST hits
            if gene_id in core_hits:
                best_hit = core_hits[gene_id]
            # If no core hit, then check the whole genome hits
            elif gene_id in genome_hits:
                best_hit = genome_hits[gene_id]
            else:
                best_hit = None
            
            # If there is a hit, write out in table
            try:
                # Parse results from best_hit without creating tuples
                hit_id = str(best_hit['sseqid'])
                coverage = str(best_hit['coverage'])
                pident = str(best_hit['pident'])
                Evalue = str(best_hit['evalue'])
                qlen = str(best_hit['qlen'])
                mismatch = str(best_hit['mismatch'])
            
                # Output into full table
                outf.write(f"{gene_id}\t{scaffold}\t{tissue_expr}\t{development_stage}\t{sex}\t{count_expr_germ}_{count_expr_soma}\t{mean_germ_TPM}_{mean_soma_TPM}\t{hit_id}\t{coverage}\t{pident}\t{Evalue}\t{qlen}\t{mismatch}\n")

            # If no hit keep original table
            except:
                outf.write(f"{gene_id}\t{scaffold}\t{tissue_expr}\t{development_stage}\t{sex}\t{count_expr_germ}_{count_expr_soma}\t{mean_germ_TPM}_{mean_soma_TPM}\tNO_BLAST_HIT\n")

# Re-order table 
def sort_table(input_file, output_file):
    # Defin order for development_stage
    stage_order = ["0-4h", "4-8h", "late-larva-early-pupa", "adult"]
    # Load the table into a pandas DataFrame
    df = pd.read_csv(input_file, sep="\t")

    # Ensure `development_stage` is treated as a categorical variable with the defined order
    df["development_stage"] = pd.Categorical(df["development_stage"], categories=stage_order, ordered=True)

    # Sort by gene_id first, then by development_stage
    df_sorted = df.sort_values(by=["gene_id", "development_stage"])

    # Save the sorted table
    df_sorted.to_csv(output_file, sep="\t", index=False)

def main():
    parser = argparse.ArgumentParser(description="Append BLAST hit info to an expression table for each GRC gene.")
    parser.add_argument("--expression_table", "-t", required=True, help="Expression table file (tab-delimited).")
    parser.add_argument("--core_blast", "-c", required=True, help="BLAST output file against the core gene list.")
    parser.add_argument("--genome_blast", "-g", required=True, help="BLAST output file against the whole genome.")
    parser.add_argument("--output", "-o", required=False, default = "GRC_BLAST_table.tsv", help="Output file for the annotated expression table.")
    args = parser.parse_args()

    get_BLAST_table(args.expression_table, args.core_blast, args.genome_blast, args.output)
    sort_table(args.output, args.output)

if __name__ == "__main__":
    main()

# python 03_get_BLAST_table.py 

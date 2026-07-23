import argparse
import os

'''
Method to define Intergenic region from https://www.biorxiv.org/content/10.1101/2022.03.31.486555v1.full
"[Intergenic regions] are defined as regions of DNA which are at least 0.5 kb from a gene annotation, 
and are at least 1 kb long. For regions which are longer than 20 kb, we keep only the 20 kb at the center 
of the region (i.e., 10 kb on each side of the center)"
'''

def read_GFF(gff_file):
    # read GFF3 file and extracts gene coordinates
    genes = []
    with open(gff_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue  # Skip comments
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue  # Skip malformed lines
            seqid, feature, start, end = fields[0], fields[2], int(fields[3]), int(fields[4])
            if feature == "gene":
                genes.append((seqid, start, end))

    # Sort genes by scaffold (seqid) and start position
    genes.sort()
    return genes

def get_intergenic_regions(genes, gff3_file, gtf_file, min_distance, min_length, max_length):
    # Identify intergenic regions and write GFF3 file
    intergenic_regions = []
    prev_end = None
    prev_seqid = None

    for seqid, start, end in genes:
        if prev_end is not None and prev_seqid == seqid:
            # Ensure minimum distance from previous gene (100bp)
            intergenic_start = prev_end + min_distance + 1
            intergenic_end = start - min_distance - 1

            # Sneure region is long enough 
            region_length = intergenic_end - intergenic_start + 1

            if region_length >= min_length:
                # Trim long regions to max_length around center
                if region_length > max_length:
                    mid_point = (intergenic_start + intergenic_end) // 2
                    intergenic_start = mid_point - (max_length // 2)
                    intergenic_end = mid_point + (max_length // 2) - 1  

                # Overlap Check: Ensure intergenic region doesn't overlap genes
                overlap = any(gene_start <= intergenic_end and gene_end >= intergenic_start 
                              for g_seqid, gene_start, gene_end in genes if g_seqid == seqid)

                if not overlap:
                    intergenic_regions.append((seqid, intergenic_start, intergenic_end))

        prev_end = end
        prev_seqid = seqid

    # Write GFF3 file
    with open(gff3_file, "w") as gff_out:
        gff_out.write("##gff-version 3\n")  # GFF3 header
        for i, (seqid, start, end) in enumerate(intergenic_regions):
            gff_out.write(f"{seqid}\tCustom\tintergenic_region\t{start}\t{end}\t.\t+\t.\tID=intergenic_{i};Name=intergenic_{i};\n")

    # Write StringTie-compatible GTF file
    with open(gtf_file, "w") as gtf_out:
        for i, (seqid, start, end) in enumerate(intergenic_regions):
            gene_id = f"intergenic_{i}"
            transcript_id = f"intergenic_{i}.t1"
            gtf_out.write(
                f"{seqid}\tCustom\texon\t{start}\t{end}\t.\t+\t.\t"
                f'gene_id "{gene_id}"; transcript_id "{transcript_id}";\n')
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate BED and GFF3 for intergenic regions from a GFF3 file")
    parser.add_argument("-g", "--gff", required=True, help="Path to annotation file (GFF3)")
    parser.add_argument("--min_distance", type=int, default=500, help="Minimum distance from genes (500 bp)")
    parser.add_argument("--min_length", type=int, default=1000, help="Minimum intergenic region length (1 kb)")
    parser.add_argument("--max_length", type=int, default=20000, help="Maximum intergenic region length (20 kb)")

    args = parser.parse_args()
    base_name = os.path.basename(args.gff).replace('.gff3', '')

    gff3_file = f"{base_name}.intergenic.gff3"
    gtf_file = f"{base_name}.intergenic.gtf"

    genes = read_GFF(args.gff)
    get_intergenic_regions(genes, gff3_file, gtf_file, args.min_distance, args.min_length, args.max_length)

# Usage: python3 01_get_intergenic_GFF3.py -g ../../Annotations/bcop_core_GRC.gff3

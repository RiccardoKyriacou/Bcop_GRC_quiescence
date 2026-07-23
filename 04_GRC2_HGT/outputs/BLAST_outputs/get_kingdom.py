import sys
import argparse
from collections import defaultdict
from Bio import Entrez

Entrez.email = "s2673271@ed.ac.uk"

def parse_all_blast_hits(blast_file):
    hits = []
    with open(blast_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            qseqid = cols[0]
            pident = cols[2]
            evalue = float(cols[10])
            staxid = cols[13]
            stitle = cols[14]

            hits.append({
                'qseqid': qseqid,
                'pident': pident,
                'evalue': evalue,
                'staxid': staxid,
                'stitle': stitle
            })
    return hits

def fetch_kingdom_for_taxids(taxids):
    """
    Queries NCBI taxonomy db for kingdoms of all taxids.
    Returns dict taxid -> kingdom
    """
    kingdoms = {}
    batch_size = 50
    taxid_list = list(taxids)
    for i in range(0, len(taxid_list), batch_size):
        batch = taxid_list[i:i+batch_size]
        ids = ",".join(batch)
        handle = Entrez.efetch(db="taxonomy", id=ids, retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        for record in records:
            taxid = record["TaxId"]
            lineage = record.get("Lineage", "")
            lineage_list = [x.strip() for x in lineage.split(";") if x.strip()]
            # Find kingdom rank: commonly 1st or 2nd lineage after "cellular organisms" or root
            kingdom = "unknown"
            for rank in lineage_list:
                if rank.lower() not in ("cellular organisms", "root"):
                    kingdom = rank
                    break
            kingdoms[taxid] = kingdom
    return kingdoms

def main():

    blast_file = "flanking_genes_BLAST_output.tsv"
    output_file = "Kingdom_info.tsv"


    all_hits = parse_all_blast_hits(blast_file)

    unique_taxids = set(hit['staxid'] for hit in all_hits)
    print(f"Fetching kingdom info for {len(unique_taxids)} unique taxids from NCBI...")
    kingdoms = fetch_kingdom_for_taxids(unique_taxids)

    print(f"Writing all hits to {output_file} ...")
    with open(output_file, "w") as out:
        for hit in all_hits:
            taxid = hit['staxid']
            kingdom = kingdoms.get(taxid, "unknown")
            gene_id = hit['qseqid']
            description = hit['stitle']
            evalue = hit['evalue']
            pident = hit['pident']
            out.write(f"{gene_id}\t{kingdom}\t{evalue}\t{pident}\t{description}\n")


if __name__ == "__main__":
    main()

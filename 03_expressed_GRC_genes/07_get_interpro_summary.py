from collections import defaultdict
import csv

input_file = "iprscan5_GRC_genes.tsv"
output_file = "interpro_GRC_summary.tsv"

annotations = defaultdict(lambda: {
    "IPR_IDs": set(),
    "IPR_Descs": set(),
    "GOs": set(),
    "Scores": [],
    "Coverages": []
})

with open(input_file, "r") as f:
    for line in f:
        if line.startswith("#") or line.strip() == "":
            continue

        fields = line.strip().split("\t")
        if len(fields) < 15:
            continue

        protein = fields[0]
        protein_length = int(fields[3]) if fields[3].isdigit() else None
        start = int(fields[6])
        end = int(fields[7])
        score = float(fields[8]) if fields[8] != "-" else None

        ipr_id = fields[11] if fields[11] != "-" else None
        ipr_desc = fields[12] if fields[12] != "-" else None
        go_terms = fields[13] if fields[13] != "-" else None
        pathways = fields[14] if fields[14] != "-" else None

        if ipr_id:
            annotations[protein]["IPR_IDs"].add(ipr_id)
        if ipr_desc:
            annotations[protein]["IPR_Descs"].add(ipr_desc)
        if go_terms:
            annotations[protein]["GOs"].update(go_terms.split("|"))

        if score is not None:
            annotations[protein]["Scores"].append(score)

        if protein_length:
            coverage = (end - start + 1) / protein_length
            annotations[protein]["Coverages"].append(coverage)

with open(output_file, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["Protein", "IPR_IDs", "IPR_Descriptions", "GO_Terms", "Best_Score", "Avg_Domain_Coverage"])

    for protein, data in annotations.items():
        best_score = min(data["Scores"]) if data["Scores"] else ""
        avg_coverage = round(sum(data["Coverages"]) / len(data["Coverages"]), 3) if data["Coverages"] else ""
        writer.writerow([
            protein,
            "; ".join(sorted(data["IPR_IDs"])),
            "; ".join(sorted(data["IPR_Descs"])),
            "; ".join(sorted(data["GOs"])),
            best_score,
            avg_coverage
        ])

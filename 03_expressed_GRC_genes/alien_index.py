#!/usr/bin/env python3
"""
alien_index.py

Computes a per-gene Alien Index (AI) from three BLASTp outputs of the same
query set against (1) a cecidomyiid proteome, (2) a GRC-excluded sciarid
core proteome, and (3) an outgroup proteome.

AI = ln(best_evalue_sciarid + 1e-200) - ln(best_evalue_cecidomyiid + 1e-200)

  AI > 0  -> better match to cecidomyiid than to sciarid
             (consistent with ancient hybridisation origin)
  AI < 0  -> better match to sciarid than to cecidomyiid
             (consistent with recent core-genome duplication)
  |AI| ~0 -> ambiguous / no strong differential signal

The outgroup best-hit is reported alongside as a sanity check: if a gene's
outgroup score is comparable to its "winning" score, that's a flag that the
AI signal may just reflect generic sequence conservation rather than a true
phylogenetic affinity.

Usage:
    python alien_index.py \
        --cecidomyiid GRC_vs_Cecidomyiid.tsv \
        --sciarid     GRC_vs_SciaridCore.tsv \
        --outgroup    GRC_vs_Outgroup.tsv \
        --out         Alien_Index_summary.csv

Expects BLAST outfmt:
    6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs
(matches the -outfmt string used in 04_run_blastp_three_way.sh)
"""

import argparse
import math
import pandas as pd

COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "qlen", "slen", "qcovs"]

NO_HIT_EVALUE = 10.0   # assigned to genes with literally zero hits in a DB
PSEUDOCOUNT = 1e-200   # avoids log(0) for very strong hits (evalue == 0.0)


def load_best_hits(path, label):
    """Read a BLAST tsv and return one row per query: its single best hit
    (highest bitscore; ties broken by lowest e-value)."""
    df = pd.read_csv(path, sep="\t", names=COLS)
    if df.empty:
        return pd.DataFrame(columns=["qseqid", f"best_evalue_{label}",
                                      f"best_bitscore_{label}",
                                      f"best_hit_{label}",
                                      f"best_pident_{label}"])
    df = df.sort_values(["qseqid", "bitscore", "evalue"],
                         ascending=[True, False, True])
    best = df.groupby("qseqid", as_index=False).first()
    best = best.rename(columns={
        "evalue": f"best_evalue_{label}",
        "bitscore": f"best_bitscore_{label}",
        "sseqid": f"best_hit_{label}",
        "pident": f"best_pident_{label}",
    })
    return best[["qseqid", f"best_evalue_{label}", f"best_bitscore_{label}",
                 f"best_hit_{label}", f"best_pident_{label}"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cecidomyiid", required=True)
    ap.add_argument("--sciarid", required=True)
    ap.add_argument("--outgroup", required=True)
    ap.add_argument("--query_fasta", required=False,
                     help="Optional: GRC_proteins.fasta, used to make sure "
                          "genes with NO hit in any DB still appear in the "
                          "output table.")
    ap.add_argument("--out", default="Alien_Index_summary.csv")
    args = ap.parse_args()

    cec = load_best_hits(args.cecidomyiid, "cecidomyiid")
    sci = load_best_hits(args.sciarid, "sciarid")
    out = load_best_hits(args.outgroup, "outgroup")

    merged = cec.merge(sci, on="qseqid", how="outer") \
                .merge(out, on="qseqid", how="outer")

    # If a query fasta was given, make sure genes with zero hits anywhere
    # still show up as rows (rather than silently disappearing).
    if args.query_fasta:
        ids = []
        with open(args.query_fasta) as f:
            for line in f:
                if line.startswith(">"):
                    ids.append(line[1:].split()[0])
        all_ids = pd.DataFrame({"qseqid": ids})
        merged = all_ids.merge(merged, on="qseqid", how="left")

    # Fill missing (no-hit) values
    for label in ["cecidomyiid", "sciarid", "outgroup"]:
        merged[f"best_evalue_{label}"] = merged[f"best_evalue_{label}"].fillna(NO_HIT_EVALUE)
        merged[f"best_bitscore_{label}"] = merged[f"best_bitscore_{label}"].fillna(0)

    # Alien Index (e-value based)
    merged["alien_index"] = merged.apply(
        lambda r: math.log(r["best_evalue_sciarid"] + PSEUDOCOUNT)
                  - math.log(r["best_evalue_cecidomyiid"] + PSEUDOCOUNT),
        axis=1
    )

    # Simple bitscore-difference version as a cross-check (less sensitive
    # to e-value saturation at very high similarity)
    merged["bitscore_diff_cec_minus_sci"] = (
        merged["best_bitscore_cecidomyiid"] - merged["best_bitscore_sciarid"]
    )

    # Flag genes where the outgroup score rivals the winning ingroup score
    # (possible generic-conservation confound rather than true affinity)
    def outgroup_flag(r):
        winner_bitscore = max(r["best_bitscore_cecidomyiid"], r["best_bitscore_sciarid"])
        if winner_bitscore == 0:
            return "no_hits_anywhere"
        ratio = r["best_bitscore_outgroup"] / winner_bitscore
        return "check_conservation" if ratio > 0.8 else "ok"

    merged["outgroup_flag"] = merged.apply(outgroup_flag, axis=1)

    merged["interpretation"] = merged["alien_index"].apply(
        lambda ai: "cecidomyiid-like" if ai > 10 else
                   ("sciarid-like" if ai < -10 else "ambiguous")
    )

    merged = merged.sort_values("alien_index", ascending=False)
    merged.to_csv(args.out, index=False)
    print(f"Wrote {args.out} ({len(merged)} genes)")
    print(merged[["qseqid", "alien_index", "interpretation", "outgroup_flag"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
import pandas as pd, re

COLS = ["qseqid","sseqid","pident","length","mismatch","gapopen",
        "qstart","qend","sstart","send","evalue","bitscore","qlen","slen"]

FILES = {
    "Bimp": "Bimp_Bcop_GRC_tblastn.tsv",
    "Ling": "Ling_Bcop_GRC_tblastn.tsv",
}

MAX_INTRON = 100_000
EVAL_CUT   = 1e-10
MIN_PID    = 60
QCOV_LEVELS = [0.7]        # produces one matrix per level


# ---------- helpers ----------
def compartment(s):
    if "GRC" in s:                 return "GRC"         
    if "GRC1" in s:                 return "GRC1"      
    if "GRC2" in s:                 return "GRC2"      
    if s.endswith("_X"):           return "X"
    if re.search(r"SUPER_\d", s):  return "core_autosome"

def merge_qcov(intervals):
    iv = sorted(intervals); out = [list(iv[0])]
    for lo, hi in iv[1:]:
        if lo <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return sum(h - l + 1 for l, h in out)

def build_loci(df):
    """Cluster alignemnts into loci per (gene, subject scaffold); return locus-level table."""
    rows = []
    for (gene, sseqid), g in df.groupby(["gene", "sseqid"]):
        g = g.copy()
        g["slo"] = g[["sstart", "send"]].min(axis=1)
        g["shi"] = g[["sstart", "send"]].max(axis=1)
        g = g.sort_values("slo")
        lid, ids, prev = 0, [], None
        for _, r in g.iterrows():
            if prev is not None and r.slo - prev > MAX_INTRON:
                lid += 1
            ids.append(lid)
            prev = r.shi if prev is None else max(prev, r.shi)
        g["locus"] = ids
        for _, sub in g.groupby("locus"):
            cov = merge_qcov(list(zip(sub.qstart, sub.qend)))
            rows.append(dict(
                gene=gene, sseqid=sseqid,
                compartment=compartment(sseqid),
                qcov=cov / sub.qlen.iloc[0],
                best_pid=sub.pident.max(),
                best_eval=sub.evalue.min(),
            ))
    return pd.DataFrame(rows)


# ---------- per-file processing ----------
def process(path):
    df = pd.read_csv(path, sep="\t", names=COLS)

    # collapse isoforms -> longest per gene
    df["gene"] = df.qseqid.str.replace(r"\.t\d+$", "", regex=True)
    keep = df.sort_values("qlen").groupby("gene").tail(1).qseqid.unique()
    df = df[df.qseqid.isin(keep)]

    # e-value gate
    df = df[df.evalue <= EVAL_CUT]

    return build_loci(df)


for sp, path in FILES.items():
    loci = process(path)

    for qcov in QCOV_LEVELS:
        good = loci[(loci.qcov >= qcov) & (loci.best_pid >= MIN_PID)]

        # count of passing loci per compartment
        n_matrix = (good.groupby(["gene", "compartment"]).size()
                        .unstack(fill_value=0))
        # best coverage per compartment (rounded for readability)
        cov_matrix = (good.groupby(["gene", "compartment"]).qcov.max()
                          .unstack().round(2))

        tag = f"{sp}_qcov{int(qcov*100)}_pident{int(MIN_PID)}"
        print(f"\n===== {sp}  (qcov>={qcov})  (pident>={MIN_PID})  — locus counts =====")
        print(n_matrix)

        n_matrix.to_csv(f"{tag}_ncounts.tsv", sep="\t")
        cov_matrix.to_csv(f"{tag}_bestqcov.tsv", sep="\t")

    # also dump the raw locus table for this species
    loci.to_csv(f"{sp}_loci_all.tsv", sep="\t", index=False)

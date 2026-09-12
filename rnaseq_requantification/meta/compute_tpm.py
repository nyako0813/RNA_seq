#!/usr/bin/env python3
"""Merge featureCounts output with gene metadata, compute TPM per sample,
and write one combined tidy CSV (+ a wide Excel workbook) comparing all
samples.

Method A (featureCounts / Subread) is used here per the instructions doc's
preference; see run_featurecounts.sh for the exact invocation
(-p --countReadPairs -t gene -g locus_tag).
"""
import pandas as pd

COUNTS_DIR = "/home/nyako/rnaseq_requantification/counts"
COUNTS_FILE = f"{COUNTS_DIR}/counts_locus_tag.txt"
GENES_FILE = "/home/nyako/rnaseq_requantification/meta/genes_table.csv"
OUT_LONG = f"{COUNTS_DIR}/rnaseq_requantification_combined.csv"
OUT_WIDE_XLSX = f"{COUNTS_DIR}/rnaseq_requantification_combined.xlsx"

# SRR -> metadata resolved from ENA/NCBI in step 1 (see meta/*.tsv and
# meta/SAMN*_biosample.xml for the underlying source records)
SAMPLE_META = {
    "SRR20281608": dict(srx="SRX16315155", biosample="SAMN29793683",
                         bioproject="PRJNA859665", substrate="acetate",
                         cell_type="humus-respiratory", replicate=3),
    "SRR20281609": dict(srx="SRX16315154", biosample="SAMN29793683",
                         bioproject="PRJNA859665", substrate="acetate",
                         cell_type="humus-respiratory", replicate=2),
    "SRR20281610": dict(srx="SRX16315153", biosample="SAMN29793683",
                         bioproject="PRJNA859665", substrate="acetate",
                         cell_type="humus-respiratory", replicate=1),
    "SRR20281611": dict(srx="SRX16315152", biosample="SAMN29793683",
                         bioproject="PRJNA859665", substrate="acetate",
                         cell_type="non-respiratory", replicate=3),
    "SRR20280855": dict(srx="SRX16314429", biosample="SAMN29787779",
                         bioproject="PRJNA859536", substrate="methanol",
                         cell_type="humus-respiratory", replicate=1),
    "SRR20280860": dict(srx="SRX16314424", biosample="SAMN29787779",
                         bioproject="PRJNA859536", substrate="methanol",
                         cell_type="humus-respiratory", replicate=2),
}

# --- load featureCounts output ---
fc = pd.read_csv(COUNTS_FILE, sep="\t", comment="#")
# columns: Geneid, Chr, Start, End, Strand, Length, <bam1>, <bam2>, ...
meta_cols = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
bam_cols = [c for c in fc.columns if c not in meta_cols]
print("BAM columns found:", bam_cols)

# map bam column -> SRR id (basename before .sorted.bam)
def bam_to_srr(col):
    base = col.split("/")[-1]
    return base.replace(".sorted.bam", "")

srr_for_col = {c: bam_to_srr(c) for c in bam_cols}
missing = [srr for srr in srr_for_col.values() if srr not in SAMPLE_META]
if missing:
    raise SystemExit(f"No metadata for SRR(s): {missing}")

genes = pd.read_csv(GENES_FILE, dtype={"old_locus_tag": str})
genes = genes.rename(columns={"locus_tag": "locus_tag"})

records = []
wide_frames = []
for col, srr in srr_for_col.items():
    m = SAMPLE_META[srr]
    sample_label = f"{srr}_{m['substrate']}_{m['cell_type']}_rep{m['replicate']}"
    sub = fc[["Geneid", "Length", col]].rename(
        columns={"Geneid": "locus_tag", col: "read_count"})
    merged = genes.merge(sub, on="locus_tag", how="left")
    merged["read_count"] = merged["read_count"].fillna(0).astype(int)
    # TPM: RPK = count / (length_kb); TPM = RPK / sum(RPK) * 1e6
    # use featureCounts' own "Length" (union of exon length for the gene
    # meta-feature) when available, else fall back to end-start+1
    length_kb = merged["Length"].fillna(
        merged["end"] - merged["start"] + 1) / 1000.0
    rpk = merged["read_count"] / length_kb
    tpm = rpk / rpk.sum() * 1e6
    merged["TPM"] = tpm
    merged["sample"] = sample_label
    merged["srr"] = srr
    for k, v in m.items():
        merged[k] = v
    records.append(merged[[
        "locus_tag", "old_locus_tag", "old_locus_tag_all", "gene_name",
        "gene_biotype", "seqid", "start", "end", "strand",
        "sample", "srr", "srx", "biosample", "bioproject", "substrate",
        "cell_type", "replicate", "read_count", "TPM",
    ]])
    w = merged[["locus_tag", "old_locus_tag", "gene_name", "gene_biotype",
                "start", "end", "strand"]].copy()
    w[f"read_count.{sample_label}"] = merged["read_count"]
    w[f"TPM.{sample_label}"] = merged["TPM"]
    wide_frames.append(w.set_index(
        ["locus_tag", "old_locus_tag", "gene_name", "gene_biotype",
         "start", "end", "strand"]))

long_df = pd.concat(records, ignore_index=True)
long_df.to_csv(OUT_LONG, index=False)
print(f"Wrote long/tidy combined CSV: {OUT_LONG} ({len(long_df)} rows)")

wide_df = pd.concat(wide_frames, axis=1).reset_index()
with pd.ExcelWriter(OUT_WIDE_XLSX, engine="openpyxl") as xw:
    wide_df.to_excel(xw, sheet_name="all_samples_wide", index=False)
    long_df.to_excel(xw, sheet_name="all_samples_long", index=False)
print(f"Wrote wide-format Excel: {OUT_WIDE_XLSX}")

# --- quick validation: known very-highly-expressed genes ---
print("\n=== Validation: rRNA / ffs / rnpB read_count & TPM per sample ===")
val = long_df[long_df["gene_biotype"].isin(["rRNA", "SRP_RNA", "RNase_P_RNA"])]
piv_counts = val.pivot_table(index=["locus_tag", "gene_name", "gene_biotype"],
                              columns="sample", values="read_count")
piv_tpm = val.pivot_table(index=["locus_tag", "gene_name", "gene_biotype"],
                           columns="sample", values="TPM")
pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 20)
print("--- read_count ---")
print(piv_counts)
print("--- TPM ---")
print(piv_tpm)

print("\n=== Max read_count / max TPM per sample (sanity check for ceiling) ===")
print(long_df.groupby("sample")[["read_count", "TPM"]].max())
print("\n=== Top 10 genes by TPM (sample 1) ===")
s1 = long_df["sample"].unique()[0]
print(long_df[long_df["sample"] == s1].sort_values("TPM", ascending=False)
      [["locus_tag", "gene_name", "gene_biotype", "read_count", "TPM"]].head(10))

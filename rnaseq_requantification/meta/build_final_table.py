#!/usr/bin/env python3
"""Merge featureCounts output with GFF3 gene metadata, compute TPM, and
write one combined CSV covering all 6 samples."""
from __future__ import annotations

import re
import urllib.parse
from pathlib import Path

import pandas as pd

GFF = Path(
    "/home/nyako/projects/ProteinHunter/data/databases/target/"
    "methanosarcina_acetivorans/ncbi_dataset/data/GCF_000007345.1/genomic.gff"
)
COUNTS = Path("/home/nyako/rnaseq_requantification/counts/counts.txt")
OUT = Path("/home/nyako/rnaseq_requantification/counts/rnaseq_requant_combined.csv")

SAMPLE_CONDITION = {
    "SRR20281608": "acetate_humus-respiratory",
    "SRR20281609": "acetate_humus-respiratory",
    "SRR20281610": "acetate_humus-respiratory",
    "SRR20281611": "acetate_no-respiratory",
    "SRR20280855": "methanol_humus-respiratory",
    "SRR20280860": "methanol_humus-respiratory",
}


def _parse_attrs(field: str) -> dict[str, str]:
    attrs = {}
    for part in field.split(";"):
        if "=" not in part:
            continue
        k, _, v = part.partition("=")
        attrs[k] = v
    return attrs


def load_gene_metadata(gff_path: Path) -> pd.DataFrame:
    rows = []
    with open(gff_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = _parse_attrs(fields[8])
            locus_tag = attrs.get("locus_tag")
            if not locus_tag:
                continue
            old_locus_tag_raw = attrs.get("old_locus_tag")
            old_locus_tag = ""
            if old_locus_tag_raw:
                decoded = urllib.parse.unquote(old_locus_tag_raw)
                parts = [p.strip() for p in decoded.split(",") if p.strip()]
                underscored = [p for p in parts if "_" in p]
                old_locus_tag = underscored[0] if underscored else (parts[0] if parts else "")
            rows.append(
                {
                    "locus_tag": locus_tag,
                    "old_locus_tag": old_locus_tag,
                    "gene_name": attrs.get("gene", ""),
                    "gene_biotype": attrs.get("gene_biotype", ""),
                    "product": attrs.get("product", ""),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    counts = pd.read_csv(COUNTS, sep="\t", comment="#")
    counts = counts.rename(columns={"Geneid": "locus_tag"})
    sample_cols = [c for c in counts.columns if c.startswith("bam/")]
    rename_map = {c: Path(c).stem.replace(".sorted", "") for c in sample_cols}
    counts = counts.rename(columns=rename_map)
    sample_names = list(rename_map.values())

    meta = load_gene_metadata(GFF)
    merged = counts.merge(meta, on="locus_tag", how="left")

    gene_length_kb = merged["Length"] / 1000.0
    for sample in sample_names:
        rpk = merged[sample] / gene_length_kb
        tpm = rpk / rpk.sum() * 1e6
        merged[f"{sample}_TPM"] = tpm
        merged = merged.rename(columns={sample: f"{sample}_read_count"})

    ordered_cols = ["locus_tag", "old_locus_tag", "gene_name", "gene_biotype", "product", "Chr", "Start", "End", "Strand", "Length"]
    for sample in sample_names:
        ordered_cols += [f"{sample}_read_count", f"{sample}_TPM"]
    merged = merged[ordered_cols]
    merged = merged.rename(columns={"Start": "start", "End": "end", "Strand": "strand", "Chr": "chr", "Length": "length_bp"})

    OUT.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(OUT, index=False)
    print(f"Wrote {OUT} ({len(merged)} genes x {len(sample_names)} samples)")

    # Sample -> condition mapping, written alongside for reference.
    cond_df = pd.DataFrame(
        [{"sample": s, "condition": SAMPLE_CONDITION.get(s, "unknown")} for s in sample_names]
    )
    cond_path = OUT.parent / "sample_conditions.csv"
    cond_df.to_csv(cond_path, index=False)
    print(f"Wrote {cond_path}")

    # Depth-cap-artifact validation: known very highly expressed genes.
    print("\n=== High-expression gene check (rRNA / ffs / rnpB) ===")
    watch_mask = (
        merged["gene_biotype"].isin(["rRNA"])
        | merged["gene_name"].isin(["ffs", "rnpB"])
        | merged["product"].str.contains("4.5S|RNase P", case=False, na=False)
    )
    watch = merged[watch_mask]
    tpm_cols = [f"{s}_TPM" for s in sample_names]
    count_cols = [f"{s}_read_count" for s in sample_names]
    print(watch[["locus_tag", "old_locus_tag", "gene_name", "gene_biotype", "product"] + count_cols].to_string(index=False))
    print("\nTPM:")
    print(watch[["locus_tag", "gene_name", "gene_biotype"] + tpm_cols].to_string(index=False))

    print("\n=== Top 15 genes by TPM (sample SRR20281608) ===")
    top = merged.sort_values("SRR20281608_TPM", ascending=False).head(15)
    print(top[["locus_tag", "old_locus_tag", "gene_name", "gene_biotype", "product", "SRR20281608_read_count", "SRR20281608_TPM"]].to_string(index=False))

    print("\n=== Max read_count per sample (sanity: is anything suspiciously capped?) ===")
    for s in sample_names:
        col = f"{s}_read_count"
        top3 = merged.nlargest(3, col)[["locus_tag", "gene_name", col]]
        print(f"-- {s} --")
        print(top3.to_string(index=False))


if __name__ == "__main__":
    main()

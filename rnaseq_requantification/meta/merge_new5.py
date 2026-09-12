#!/usr/bin/env python3
"""Merge featureCounts output for the 5 additional SRA runs into the
existing rnaseq_requant_combined.xlsx/.csv, using the same TPM formula as
meta/build_final_table.py (RPK = count / (Length_kb); TPM = RPK / sum(RPK) * 1e6,
using featureCounts' own gene-level "Length" column)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path("/home/nyako/rnaseq_requantification")
COUNTS_NEW = ROOT / "counts/counts_new5.txt"
XLSX = ROOT / "counts/rnaseq_requant_combined.xlsx"
CSV = ROOT / "counts/rnaseq_requant_combined.csv"

NEW_SAMPLE_CONDITION = {
    "SRR20281612": "acetate_no-respiratory",
    "SRR20281613": "acetate_no-respiratory",
    "SRR20650029": "methanol_no-respiratory",
    "SRR20650030": "methanol_no-respiratory",
    "SRR20650031": "methanol_no-respiratory",
}


def main() -> None:
    new_counts = pd.read_csv(COUNTS_NEW, sep="\t", comment="#")
    new_counts = new_counts.rename(columns={"Geneid": "locus_tag"})
    bam_cols = [c for c in new_counts.columns if c.startswith("bam/")]
    rename_map = {c: Path(c).stem.replace(".sorted", "") for c in bam_cols}
    new_counts = new_counts.rename(columns=rename_map)
    new_samples = list(rename_map.values())
    print("New samples:", new_samples)

    gene_length_kb = new_counts["Length"] / 1000.0
    for sample in new_samples:
        rpk = new_counts[sample] / gene_length_kb
        tpm = rpk / rpk.sum() * 1e6
        new_counts[f"{sample}_TPM"] = tpm
        new_counts = new_counts.rename(columns={sample: f"{sample}_read_count"})

    new_cols = ["locus_tag"]
    for sample in new_samples:
        new_cols += [f"{sample}_read_count", f"{sample}_TPM"]
    new_slim = new_counts[new_cols]

    # --- existing gene_counts_tpm sheet ---
    existing = pd.read_excel(XLSX, sheet_name="gene_counts_tpm")
    n_before = len(existing)
    merged = existing.merge(new_slim, on="locus_tag", how="left", validate="one_to_one")
    assert len(merged) == n_before, f"row count changed: {n_before} -> {len(merged)}"
    missing = merged[[f"{s}_read_count" for s in new_samples]].isna().any(axis=1)
    assert not missing.any(), f"{missing.sum()} genes failed to match by locus_tag"
    for s in new_samples:
        merged[f"{s}_read_count"] = merged[f"{s}_read_count"].astype(int)

    # --- existing sample_conditions sheet ---
    existing_cond = pd.read_excel(XLSX, sheet_name="sample_conditions")
    new_cond_rows = pd.DataFrame(
        [{"sample": s, "condition": NEW_SAMPLE_CONDITION[s]} for s in new_samples]
    )
    merged_cond = pd.concat([existing_cond, new_cond_rows], ignore_index=True)

    with pd.ExcelWriter(XLSX, engine="openpyxl") as xw:
        merged.to_excel(xw, sheet_name="gene_counts_tpm", index=False)
        merged_cond.to_excel(xw, sheet_name="sample_conditions", index=False)
    print(f"Updated {XLSX}: {len(merged)} genes x "
          f"{(len(merged.columns) - 10) // 2} samples")

    merged.to_csv(CSV, index=False)
    merged_cond.to_csv(CSV.parent / "sample_conditions.csv", index=False)
    print(f"Updated {CSV} and sample_conditions.csv")

    # --- validation ---
    print("\n=== TPM column sums (should be ~1e6) ===")
    for s in new_samples:
        print(f"  {s}: {merged[f'{s}_TPM'].sum():,.2f}")

    print("\n=== High-expression gene check (rRNA / ffs / rnpB) for new samples ===")
    watch_mask = (
        merged["gene_biotype"].isin(["rRNA"])
        | merged["gene_name"].isin(["ffs", "rnpB"])
    )
    watch = merged[watch_mask]
    cols = ["locus_tag", "gene_name", "gene_biotype"]
    for s in new_samples:
        cols += [f"{s}_read_count", f"{s}_TPM"]
    print(watch[cols].to_string(index=False))


if __name__ == "__main__":
    main()

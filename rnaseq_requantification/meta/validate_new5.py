#!/usr/bin/env python3
"""Step 4 validation: mapping-rate summary, TPM-sum sanity check (done in
merge_new5.py already), and within-condition replicate correlation
(Pearson on log2(TPM+1)) for the two now-complete triplicates."""
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/nyako/rnaseq_requantification")
BAM_DIR = ROOT / "bam"
CSV = ROOT / "counts/rnaseq_requant_combined.csv"

ALL_SAMPLES = [
    "SRR20281608", "SRR20281609", "SRR20281610", "SRR20281611",
    "SRR20280855", "SRR20280860",
    "SRR20281612", "SRR20281613", "SRR20650029", "SRR20650030", "SRR20650031",
]

print("=== Mapping-rate summary (from samtools flagstat) ===")
rows = []
for s in ALL_SAMPLES:
    fs = BAM_DIR / f"{s}.flagstat.txt"
    if not fs.exists():
        rows.append({"sample": s, "total": None, "mapped_pct": None, "properly_paired_pct": None})
        continue
    text = fs.read_text()
    total_m = re.search(r"^(\d+) \+ \d+ in total", text, re.M)
    mapped_m = re.search(r"^\d+ \+ \d+ mapped \(([\d.]+)%", text, re.M)
    pp_m = re.search(r"^\d+ \+ \d+ properly paired \(([\d.]+)%", text, re.M)
    rows.append({
        "sample": s,
        "total_alignments": int(total_m.group(1)) if total_m else None,
        "mapped_pct": float(mapped_m.group(1)) if mapped_m else None,
        "properly_paired_pct": float(pp_m.group(1)) if pp_m else None,
    })
summary = pd.DataFrame(rows)
print(summary.to_string(index=False))

print("\n=== Within-condition replicate correlation (Pearson, log2(TPM+1)) ===")
df = pd.read_csv(CSV)

groups = {
    "acetate_no-respiratory (611/612/613)": ["SRR20281611", "SRR20281612", "SRR20281613"],
    "methanol_no-respiratory (029/030/031)": ["SRR20650029", "SRR20650030", "SRR20650031"],
}
for label, samples in groups.items():
    print(f"\n-- {label} --")
    log_tpm = pd.DataFrame({s: np.log2(df[f"{s}_TPM"] + 1) for s in samples})
    corr = log_tpm.corr(method="pearson")
    print(corr.round(4).to_string())

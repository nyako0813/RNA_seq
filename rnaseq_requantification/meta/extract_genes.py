#!/usr/bin/env python3
"""Extract a gene coordinate/attribute table from the M. acetivorans GFF3.

Read-only source (outside ProteinHunter repo policy: this script only READS
the GFF, never modifies it):
/home/nyako/projects/ProteinHunter/data/databases/target/methanosarcina_acetivorans/ncbi_dataset/data/GCF_000007345.1/genomic.gff
"""
import re
import urllib.parse
import pandas as pd

GFF = "/home/nyako/projects/ProteinHunter/data/databases/target/methanosarcina_acetivorans/ncbi_dataset/data/GCF_000007345.1/genomic.gff"
OUT = "/home/nyako/rnaseq_requantification/meta/genes_table.csv"


def parse_attrs(field):
    d = {}
    for kv in field.strip().split(";"):
        if not kv:
            continue
        if "=" not in kv:
            continue
        k, v = kv.split("=", 1)
        d[k] = urllib.parse.unquote(v)
    return d


rows = []
with open(GFF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        seqid, source, ftype, start, end, score, strand, frame, attrs = parts
        if ftype != "gene":
            continue
        a = parse_attrs(attrs)
        locus_tag = a.get("locus_tag", "")
        old_locus_tag_raw = a.get("old_locus_tag", "")
        # old_locus_tag can be a comma-separated list e.g. "MA0001,MA_0001"
        # keep the raw (comma-joined) value AND a normalized "MA_####"-style
        # first token for convenience, matching the convention used
        # elsewhere in the adjacent ProteinHunter Rockhopper work.
        old_locus_tag_list = [x for x in old_locus_tag_raw.split(",") if x]
        old_locus_tag_primary = old_locus_tag_list[0] if old_locus_tag_list else ""
        gene_name = a.get("gene", a.get("Name", ""))
        gene_biotype = a.get("gene_biotype", "")
        rows.append({
            "locus_tag": locus_tag,
            "old_locus_tag": old_locus_tag_primary,
            "old_locus_tag_all": old_locus_tag_raw,
            "gene_name": gene_name,
            "gene_biotype": gene_biotype,
            "seqid": seqid,
            "start": int(start),
            "end": int(end),
            "strand": strand,
        })

df = pd.DataFrame(rows)
df.to_csv(OUT, index=False)
print(f"Wrote {len(df)} gene records to {OUT}")
print(df["gene_biotype"].value_counts())
print("Missing old_locus_tag:", (df["old_locus_tag"] == "").sum())

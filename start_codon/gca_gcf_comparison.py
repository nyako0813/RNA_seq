#!/usr/bin/env python3
"""
gca_gcf_comparison.py

目的
----
現行パイプラインが使っているRefSeqアノテーション(GCF_000007345.1、
2024年10月版)と、2002年の原論文由来のGenBankアノテーション
(GCA_000007345.1)を比較し、両者で翻訳開始点(開始コドン位置)が
異なる遺伝子を洗い出す。

ゲノム配列自体はAE010299.1(GCA)とNC_003552.1(GCF)で完全一致することを
別途確認済み(diffなし、5,751,492bp)。従って座標はそのまま比較できる。

マッチング方法
--------------
GCFの`gene`行の`old_locus_tag`属性(例: "MA0001,MA_0001")から
"MA_XXXX"形式を取り出し、GCAの`gene`行の`locus_tag`(そのままMA_XXXX形式)
と突き合わせる。両方とも`gene_biotype=protein_coding`の行のみを対象とする
(tRNA/rRNA/pseudogeneには開始コドンの概念がなじまないため)。

使い方
------
    python gca_gcf_comparison.py \\
        --gcf-gff <GCF gff> \\
        --gca-gff <GCA gff> \\
        --out full_run/gca_gcf_comparison.csv
"""

from __future__ import annotations

import argparse
import csv
import re
import urllib.parse
from dataclasses import dataclass
from typing import Dict, List, Optional


@dataclass
class GeneRecord:
    locus_tag: str
    start: int
    end: int
    strand: str
    old_locus_tag: Optional[str] = None


def parse_attrs(attr_str: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for kv in attr_str.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            attrs[k] = v
    return attrs


def parse_old_locus_tag(raw: Optional[str]) -> Optional[str]:
    """"MA0001%2CMA_0001"のようなold_locus_tag属性から"MA_XXXX"形式を取り出す。

    find_true_start_codons.py の同名関数とロジックを揃えている。
    """
    if not raw:
        return None
    decoded = urllib.parse.unquote(raw)
    candidates = [v.strip() for v in decoded.split(",") if v.strip()]
    for c in candidates:
        if re.fullmatch(r"MA_\d+", c):
            return c
    return decoded


def parse_protein_coding_genes(path: str) -> Dict[str, GeneRecord]:
    """GFFの`gene`行のうちgene_biotype=protein_codingのものだけをlocus_tagで引けるようにする。"""
    genes: Dict[str, GeneRecord] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9 or cols[2] != "gene":
                continue
            seqid, source, ftype, start, end, score, strand, frame, attr_str = cols
            attrs = parse_attrs(attr_str)
            if attrs.get("gene_biotype") != "protein_coding":
                continue
            locus_tag = attrs.get("locus_tag")
            if not locus_tag:
                continue
            genes[locus_tag] = GeneRecord(
                locus_tag=locus_tag,
                start=int(start),
                end=int(end),
                strand=strand,
                old_locus_tag=parse_old_locus_tag(attrs.get("old_locus_tag")),
            )
    return genes


def translation_start(gene: GeneRecord) -> int:
    """+鎖ならstart、-鎖ならendを翻訳開始位置として返す。"""
    return gene.start if gene.strand == "+" else gene.end


def build_comparison(gcf_genes: Dict[str, GeneRecord], gca_genes: Dict[str, GeneRecord]) -> List[dict]:
    # GCAはlocus_tagがそのままMA_XXXX形式なのでold_locus_tagとの突合キーになる
    gca_by_locus = gca_genes

    rows: List[dict] = []
    for locus_tag, gcf_gene in sorted(gcf_genes.items()):
        old_locus_tag = gcf_gene.old_locus_tag
        if not old_locus_tag or old_locus_tag not in gca_by_locus:
            continue
        gca_gene = gca_by_locus[old_locus_tag]
        if gcf_gene.strand != gca_gene.strand:
            # strandが食い違う場合は座標比較の前提が崩れるため個別記録のみ行い、diffは計算しない
            rows.append({
                "locus_tag": locus_tag,
                "old_locus_tag": old_locus_tag,
                "gcf_start": translation_start(gcf_gene),
                "gca_start": translation_start(gca_gene),
                "diff_bp": "",
                "diff_codons": "",
                "frame_shift_flag": "",
                "strand_mismatch": True,
            })
            continue

        gcf_start = translation_start(gcf_gene)
        gca_start = translation_start(gca_gene)

        if gcf_gene.strand == "+":
            # +鎖: 座標が小さいほど上流(N末端側)
            diff_bp = gcf_start - gca_start
        else:
            # -鎖: 座標が大きいほど上流(N末端側)
            diff_bp = gca_start - gcf_start
        # diff_bp > 0 : GCAが現行(GCF)よりも上流の開始点を示唆(延長方向)
        # diff_bp < 0 : GCAが現行(GCF)よりも下流の開始点を示唆(短縮方向)

        frame_shift = (diff_bp % 3) != 0
        diff_codons = diff_bp / 3

        rows.append({
            "locus_tag": locus_tag,
            "old_locus_tag": old_locus_tag,
            "gcf_start": gcf_start,
            "gca_start": gca_start,
            "diff_bp": diff_bp,
            "diff_codons": diff_codons,
            "frame_shift_flag": frame_shift,
            "strand_mismatch": False,
        })
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gcf-gff", required=True)
    ap.add_argument("--gca-gff", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    gcf_genes = parse_protein_coding_genes(args.gcf_gff)
    gca_genes = parse_protein_coding_genes(args.gca_gff)
    print(f"GCF protein_coding gene数: {len(gcf_genes)}")
    print(f"GCA protein_coding gene数: {len(gca_genes)}")

    rows = build_comparison(gcf_genes, gca_genes)
    matched = [r for r in rows if not r["strand_mismatch"]]
    mismatched_strand = [r for r in rows if r["strand_mismatch"]]
    diff_rows = [r for r in matched if r["diff_bp"] != 0]
    frame_shift_rows = [r for r in matched if r["frame_shift_flag"] is True]

    print(f"old_locus_tag経由でマッチした遺伝子数: {len(rows)}")
    if mismatched_strand:
        print(f"  うちstrand不一致(座標比較スキップ): {len(mismatched_strand)}件")
    print(f"開始位置が一致: {len(matched) - len(diff_rows)}件")
    print(f"開始位置が不一致: {len(diff_rows)}件")
    print(f"  うちフレームシフト疑い(差分が3の倍数でない): {len(frame_shift_rows)}件")

    fieldnames = ["locus_tag", "old_locus_tag", "gcf_start", "gca_start", "diff_bp", "diff_codons", "frame_shift_flag"]
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r[k] for k in fieldnames})
    print(f"出力: {args.out}")


if __name__ == "__main__":
    main()

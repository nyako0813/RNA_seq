#!/usr/bin/env python3
"""
resume_gca_gcf_step9e.py

目的
----
`start_codon_analysis_summary.md` §9-E で一時停止していたGCA/GCF比較の
追加検証(Step1〜3)を実行する。PC移行後、RNA-seq BAM・オルソログ蛋白質
BLAST DB・samtoolsが揃ったため、新規のBLAST/RNA-seq計算を行う。

Step 1: フレームシフト疑い6件について、GCF読み枠・GCA読み枠それぞれで
        翻訳したタンパク質をオルソログDBにBLASTし、どちらの読み枠が
        妥当かを判定する(`gca_frameshift_genes.csv`を上書き更新)。
Step 2: `gca_flagged_341.csv`のうち`weak_evidence_exclude=False`の207件
        について、GCA提案位置(`gca_start`)での新規BLAST同一性・
        RNA-seqステップ比を計算し、既存パイプラインの判定基準
        (ステップ比閾値3.0)に沿って再分類する。
Step 3: `gca_flagged_850.csv`のうち「GCAが一致せず中確信度を維持」の
        119件について、Step2と同様の新規計算を行い再分類する。

判定基準(既存パイプラインと揃えた値)
--------------------------------
- RNA-seqステップ比の閾値: 3.0 (find_true_start_codons.COVERAGE_STEP_THRESHOLDと同じ)
- ステップ比計算窓: 30bp (find_true_start_codons.COVERAGE_WINDOWと同じ)
- BLASTの「N末端付近から良いアラインメント」とみなす基準: qstart <= 3
- BLASTヒットを「弱すぎない」とみなす同一性の下限: pident >= 25.0
- GCA提案位置の候補コドンが非正準、または現行ストップコドンとの間に
  読み枠内ストップコドンがある場合は、BLAST/RNA-seqの数値がどうであれ
  「支持」とは認めない(§9-Cで既存の341件フラグ付けに用いた
  `check_upstream_extension_valid`と同じ考え方)。

使い方
------
    python resume_gca_gcf_step9e.py
"""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

BASE = Path(__file__).resolve().parent
FULL_RUN = BASE / "full_run"
REPO_ROOT = BASE.parent

GENOME_FASTA = REPO_ROOT / "rnaseq_requantification" / "genome" / "NC_003552.1.fna"
GCF_GFF = BASE / "data" / "methanosarcina_acetivorans" / "ncbi_dataset" / "data" / "GCF_000007345.1" / "genomic.gff"
ORTHOLOGS_FAA = FULL_RUN / "orthologs_methanosarcina.faa"
SAMTOOLS_BIN = REPO_ROOT / "rnaseq_requantification" / "tools" / "mamba_root" / "envs" / "rnaseq" / "bin" / "samtools"
BAM_DIR = REPO_ROOT / "rnaseq_requantification" / "bam"
BAM_PATHS = sorted(BAM_DIR.glob("*.sorted.bam"))

SEQID = "NC_003552.1"
WORKDIR = FULL_RUN / ".step9e_work"

sys.path.insert(0, str(BASE))
from gca_gcf_comparison import parse_protein_coding_genes, GeneRecord  # noqa: E402
from integrate_gca_evidence import (  # noqa: E402
    load_fasta, revcomp, translate, get_orf_protein, START_CODONS,
)
from find_true_start_codons import (  # noqa: E402
    get_pooled_depth, mean_depth, run_blast, parse_best_blast_hits,
    COVERAGE_WINDOW, COVERAGE_STEP_THRESHOLD,
)

BLAST_QSTART_GOOD = 3
BLAST_PIDENT_WEAK = 25.0


# ---------------------------------------------------------------------------
# 共通ユーティリティ
# ---------------------------------------------------------------------------

def orf_protein_forced_met(seq: str, lo: int, hi: int, strand: str) -> str:
    """開始コドンをMとして翻訳する(find_true_start_codons.translate_cdsと同じ規約)。"""
    aa = get_orf_protein(seq, lo, hi, strand)
    if aa and aa[0] != "M":
        aa = "M" + aa[1:]
    return aa


def strip_stop(aa: str) -> str:
    return aa[:-1] if aa.endswith("*") else aa


def run_blastp_batch(query_fasta: Path, out_tsv: Path) -> Dict[str, dict]:
    if not query_fasta.exists() or query_fasta.stat().st_size == 0:
        return {}
    run_blast(str(query_fasta), str(ORTHOLOGS_FAA), str(out_tsv), threads=4)
    return parse_best_blast_hits(str(out_tsv))


def step_ratio_at(depth: Dict[int, int], pos: int, strand: str, window: int = COVERAGE_WINDOW) -> float:
    """find_true_start_codons.score_candidates_by_coverageと同じ式を単一候補用に適用する。"""
    if strand == "+":
        up_start, up_end = pos - window, pos - 1
        down_start, down_end = pos, pos + window - 1
    else:
        up_start, up_end = pos + 1, pos + window
        down_start, down_end = pos - window + 1, pos
    up_mean = mean_depth(depth, up_start, up_end)
    down_mean = mean_depth(depth, down_start, down_end)
    return (down_mean + 1.0) / (up_mean + 1.0)


def rna_step_ratio_at_pos(pos: int, strand: str) -> float:
    pad = COVERAGE_WINDOW + 5
    region_start = max(1, pos - pad)
    region_end = pos + pad
    depth = get_pooled_depth([str(p) for p in BAM_PATHS], SEQID, region_start, region_end, str(SAMTOOLS_BIN))
    return round(step_ratio_at(depth, pos, strand), 2)


# ---------------------------------------------------------------------------
# Step 1: フレームシフト疑い6件
# ---------------------------------------------------------------------------

def recommend_frame(gcf_hit: Optional[dict], gca_hit: Optional[dict]) -> str:
    def is_good(hit: Optional[dict]) -> bool:
        return hit is not None and hit["qstart"] <= BLAST_QSTART_GOOD and hit["pident"] >= BLAST_PIDENT_WEAK

    gcf_good = is_good(gcf_hit)
    gca_good = is_good(gca_hit)

    if gcf_good and not gca_good:
        return "GCF読み枠が妥当"
    if gca_good and not gcf_good:
        return "GCA読み枠が妥当"
    if gcf_good and gca_good:
        if abs(gcf_hit["pident"] - gca_hit["pident"]) < 5.0 and gcf_hit["qstart"] == gca_hit["qstart"]:
            return "判断つかない(両方ヒットが弱い、または拮抗)"
        return "GCF読み枠が妥当" if gcf_hit["pident"] >= gca_hit["pident"] else "GCA読み枠が妥当"
    return "判断つかない(両方ヒットが弱い、または拮抗)"


def step1_frameshift(seq: str) -> None:
    in_path = FULL_RUN / "gca_frameshift_genes.csv"
    with open(in_path) as fh:
        rows = list(csv.DictReader(fh))

    WORKDIR.mkdir(parents=True, exist_ok=True)
    query_fasta = WORKDIR / "step1_frame_proteins.faa"
    with open(query_fasta, "w") as fh:
        for row in rows:
            strand = row["strand"]
            gcf_start, gcf_end = int(row["gcf_start"]), int(row["gcf_end"])
            gca_start, gca_end = int(row["gca_start"]), int(row["gca_end"])
            gcf_lo, gcf_hi = (gcf_start, gcf_end) if strand == "+" else (gcf_end, gcf_start)
            gca_lo, gca_hi = (gca_start, gca_end) if strand == "+" else (gca_end, gca_start)
            gcf_prot = strip_stop(orf_protein_forced_met(seq, gcf_lo, gcf_hi, strand))
            gca_prot = strip_stop(orf_protein_forced_met(seq, gca_lo, gca_hi, strand))
            fh.write(f">{row['locus_tag']}__gcf\n{gcf_prot}\n")
            fh.write(f">{row['locus_tag']}__gca\n{gca_prot}\n")

    best = run_blastp_batch(query_fasta, WORKDIR / "step1_blast.tsv")

    out_rows = []
    for row in rows:
        locus_tag = row["locus_tag"]
        gcf_hit = best.get(f"{locus_tag}__gcf")
        gca_hit = best.get(f"{locus_tag}__gca")
        out = dict(row)
        out["gcf_frame_best_pident"] = gcf_hit["pident"] if gcf_hit else ""
        out["gcf_frame_best_qstart"] = gcf_hit["qstart"] if gcf_hit else ""
        out["gca_frame_best_pident"] = gca_hit["pident"] if gca_hit else ""
        out["gca_frame_best_qstart"] = gca_hit["qstart"] if gca_hit else ""
        out["recommended_frame"] = recommend_frame(gcf_hit, gca_hit)
        out_rows.append(out)

    fieldnames = list(rows[0].keys()) + [
        "gcf_frame_best_pident", "gcf_frame_best_qstart",
        "gca_frame_best_pident", "gca_frame_best_qstart", "recommended_frame",
    ]
    with open(in_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"[Step1] {len(out_rows)}件を{in_path}に上書き出力しました。")
    for r in out_rows:
        print(
            f"    {r['locus_tag']}: {r['recommended_frame']} "
            f"(GCF: pident={r['gcf_frame_best_pident']}, qstart={r['gcf_frame_best_qstart']} / "
            f"GCA: pident={r['gca_frame_best_pident']}, qstart={r['gca_frame_best_qstart']})"
        )


# ---------------------------------------------------------------------------
# Step 2 / Step 3 共通: GCA提案位置での新規BLAST・RNA-seq証拠
# ---------------------------------------------------------------------------

def orf_span_from_gca_start(gene: GeneRecord, gca_start: int) -> Tuple[int, int]:
    """gca_start(候補開始点)からGCF側と共通のストップコドンまでのORF区間(low, high)を返す。
    非フレームシフト前提(このスクリプトが対象とする341/850系列はいずれもフレームシフト除外済み)
    のため、ストップコドン位置はGCF/GCA間で共通。
    """
    if gene.strand == "+":
        return gca_start, gene.end
    return gene.start, gca_start


def check_gca_start_validity(seq: str, lo: int, hi: int, strand: str) -> Tuple[bool, str]:
    """候補コドンが正準(ATG/GTG/TTG)か、区間内に読み枠内ストップコドンが無いかを確認する。
    (§9-Cの`check_upstream_extension_valid`と同じ考え方を、延長・短縮どちらの方向にも適用する形に一般化)
    """
    if strand == "+":
        candidate_codon = seq[lo - 1: lo + 2]
    else:
        candidate_codon = revcomp(seq[hi - 3: hi])
    if candidate_codon not in START_CODONS:
        return False, f"候補位置のコドンが非正準({candidate_codon})"
    aa = get_orf_protein(seq, lo, hi, strand)
    if "*" in aa[:-1]:
        return False, "候補位置と現行ストップコドンの間に読み枠内ストップコドンあり"
    return True, "妥当"


def load_gcf_genes() -> Dict[str, GeneRecord]:
    return parse_protein_coding_genes(str(GCF_GFF))


def compute_new_evidence_for_rows(
    rows: List[dict], gcf_genes: Dict[str, GeneRecord], seq: str, tag: str
) -> Dict[str, dict]:
    WORKDIR.mkdir(parents=True, exist_ok=True)
    query_fasta = WORKDIR / f"{tag}_proteins.faa"

    spans: Dict[str, Tuple[int, int, str]] = {}
    validity: Dict[str, Tuple[bool, str]] = {}
    with open(query_fasta, "w") as fh:
        for row in rows:
            locus_tag = row["locus_tag"]
            gene = gcf_genes[locus_tag]
            gca_start = int(float(row["gca_start"]))
            lo, hi = orf_span_from_gca_start(gene, gca_start)
            spans[locus_tag] = (lo, hi, gene.strand)
            valid, reason = check_gca_start_validity(seq, lo, hi, gene.strand)
            validity[locus_tag] = (valid, reason)
            if valid:
                prot = strip_stop(orf_protein_forced_met(seq, lo, hi, gene.strand))
                fh.write(f">{locus_tag}\n{prot}\n")

    best = run_blastp_batch(query_fasta, WORKDIR / f"{tag}_blast.tsv")

    evidence: Dict[str, dict] = {}
    for i, row in enumerate(rows, 1):
        locus_tag = row["locus_tag"]
        lo, hi, strand = spans[locus_tag]
        gca_start = int(float(row["gca_start"]))
        valid, reason = validity[locus_tag]
        hit = best.get(locus_tag)
        blast_pident = hit["pident"] if hit else None
        blast_qstart = hit["qstart"] if hit else None
        blast_supports = (
            valid and hit is not None
            and blast_qstart <= BLAST_QSTART_GOOD and blast_pident >= BLAST_PIDENT_WEAK
        )
        rna_ratio = rna_step_ratio_at_pos(gca_start, strand)
        rna_supports = valid and rna_ratio >= COVERAGE_STEP_THRESHOLD
        evidence[locus_tag] = {
            "orf_valid_at_gca": valid,
            "orf_note": reason,
            "blast_pident_at_gca": blast_pident,
            "blast_qstart_at_gca": blast_qstart,
            "rna_step_ratio_at_gca": rna_ratio,
            "blast_supports_gca": blast_supports,
            "rna_supports_gca": rna_supports,
        }
        print(
            f"  [{tag} {i}/{len(rows)}] {locus_tag}: valid={valid}"
            f" blast_pident={blast_pident} qstart={blast_qstart} rna_ratio={rna_ratio}",
            file=sys.stderr,
        )
    return evidence


# ---------------------------------------------------------------------------
# Step 2: 「支持」336件中、根拠薄弱と判定されなかった207件
# ---------------------------------------------------------------------------

def step2(gcf_genes: Dict[str, GeneRecord], seq: str) -> List[dict]:
    in_path = FULL_RUN / "gca_flagged_341.csv"
    with open(in_path) as fh:
        rows = list(csv.DictReader(fh))
    target_rows = [r for r in rows if r["weak_evidence_exclude"] == "False"]
    print(f"[Step2] 対象: {len(target_rows)}件(gca_flagged_341.csvのweak_evidence_exclude=False)")

    evidence = compute_new_evidence_for_rows(target_rows, gcf_genes, seq, tag="step2")

    out_rows = []
    for row in target_rows:
        ev = evidence[row["locus_tag"]]
        if ev["blast_supports_gca"] and ev["rna_supports_gca"]:
            reclass = "高確信度相当に格上げ(新規BLAST・RNA-seq証拠の両方がGCA提案位置を支持)"
        elif ev["blast_supports_gca"] or ev["rna_supports_gca"]:
            reclass = "中確信度相当に格上げ(新規証拠の片方のみがGCA提案位置を支持)"
        else:
            reclass = "現行アノテーションを支持のまま(新規証拠はGCA提案位置を支持せず)"
        out = dict(row)
        out["blast_pident_at_gca"] = ev["blast_pident_at_gca"] if ev["blast_pident_at_gca"] is not None else ""
        out["blast_qstart_at_gca"] = ev["blast_qstart_at_gca"] if ev["blast_qstart_at_gca"] is not None else ""
        out["rna_step_ratio_at_gca"] = ev["rna_step_ratio_at_gca"]
        out["orf_valid_at_gca"] = ev["orf_valid_at_gca"]
        out["orf_note"] = ev["orf_note"]
        out["reclassification"] = reclass
        out_rows.append(out)

    out_path = FULL_RUN / "gca_flagged_341_reevaluated.csv"
    fieldnames = list(rows[0].keys()) + [
        "blast_pident_at_gca", "blast_qstart_at_gca", "rna_step_ratio_at_gca",
        "orf_valid_at_gca", "orf_note", "reclassification",
    ]
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(out_rows)

    counts = Counter(r["reclassification"].split("(")[0] for r in out_rows)
    print(f"[Step2] 出力: {out_path} ({len(out_rows)}件)")
    for k, v in counts.items():
        print(f"    {k}: {v}件")
    return out_rows


# ---------------------------------------------------------------------------
# Step 3: 「中確信度」849件中、GCAと不一致の119件
# ---------------------------------------------------------------------------

def step3(gcf_genes: Dict[str, GeneRecord], seq: str) -> List[dict]:
    in_path = FULL_RUN / "gca_flagged_850.csv"
    with open(in_path) as fh:
        rows = list(csv.DictReader(fh))
    target_rows = [r for r in rows if not r["revised_confidence"].startswith("高確信度相当")]
    print(f"[Step3] 対象: {len(target_rows)}件(gca_flagged_850.csvのGCAが一致せず中確信度を維持)")

    evidence = compute_new_evidence_for_rows(target_rows, gcf_genes, seq, tag="step3")

    out_rows = []
    for row in target_rows:
        ev = evidence[row["locus_tag"]]
        if ev["blast_supports_gca"] and ev["rna_supports_gca"]:
            reclass = "高確信度相当への格上げ(新規BLAST・RNA-seq証拠の両方がGCA提案位置を支持)"
        elif ev["blast_supports_gca"]:
            reclass = "高確信度相当への格上げ(新規BLAST証拠がGCA提案位置を支持)"
        elif ev["rna_supports_gca"]:
            reclass = "高確信度相当への格上げ(新規RNA-seq証拠がGCA提案位置を支持)"
        else:
            reclass = "中確信度のまま(新規証拠もGCA提案位置を支持せず)"
        out = dict(row)
        out["blast_pident_at_gca"] = ev["blast_pident_at_gca"] if ev["blast_pident_at_gca"] is not None else ""
        out["blast_qstart_at_gca"] = ev["blast_qstart_at_gca"] if ev["blast_qstart_at_gca"] is not None else ""
        out["rna_step_ratio_at_gca"] = ev["rna_step_ratio_at_gca"]
        out["orf_valid_at_gca"] = ev["orf_valid_at_gca"]
        out["orf_note"] = ev["orf_note"]
        out["reclassification"] = reclass
        out_rows.append(out)

    out_path = FULL_RUN / "gca_flagged_850_reevaluated.csv"
    fieldnames = list(rows[0].keys()) + [
        "blast_pident_at_gca", "blast_qstart_at_gca", "rna_step_ratio_at_gca",
        "orf_valid_at_gca", "orf_note", "reclassification",
    ]
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(out_rows)

    upgraded = sum(1 for r in out_rows if r["reclassification"].startswith("高確信度相当"))
    print(f"[Step3] 出力: {out_path} ({len(out_rows)}件中 {upgraded}件を高確信度相当へ格上げ提案)")
    return out_rows


# ---------------------------------------------------------------------------
# メイン処理
# ---------------------------------------------------------------------------

def main() -> None:
    if not BAM_PATHS:
        raise SystemExit(f"BAMファイルが見つかりません: {BAM_DIR}")
    print(f"BAM: {len(BAM_PATHS)}件 ({', '.join(p.name for p in BAM_PATHS)})")

    genome = load_fasta(GENOME_FASTA)
    seq = genome[SEQID]
    print(f"ゲノム配列読み込み: {SEQID} ({len(seq)}bp)")

    step1_frameshift(seq)

    gcf_genes = load_gcf_genes()
    print(f"GCF protein_coding gene数: {len(gcf_genes)}")

    step2(gcf_genes, seq)
    step3(gcf_genes, seq)


if __name__ == "__main__":
    main()

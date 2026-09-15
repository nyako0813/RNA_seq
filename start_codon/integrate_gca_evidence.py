#!/usr/bin/env python3
"""
integrate_gca_evidence.py

目的
----
`full_results.csv`(RNA-seq・BLASTによる開始コドン再検証パイプラインの
全4,683遺伝子の結果)に、2002年原論文由来のGenBankアノテーション(GCA)の
開始位置(`gca_start`)を突き合わせ、以下を行う。

  1. `gca_start`が`rna_best_pos`・`blast_supported_pos`・現行開始位置の
     どれと一致するか(`gca_agrees_with`)を判定する。
  2. その一致状況をもとに、既存の`verdict`(および「要確認」34件については
     `review_34_final_merged.csv`の`final_classification`)を出発点として
     `revised_confidence`を提案する。

注意(実行環境の制約)
--------------------
本来のStep 3-B/3-Cは、GCAが提案する新規候補位置についてBLAST/RNA-seqを
「新たに計算」することを求めているが、当時使用したBAMファイル・オルソログ
DB・samtoolsは本環境に現存しないため、**新規のBLAST/RNA-seq計算は行わない**
(ユーザー確認済み)。代わりに、既存の`full_results.csv`の列
(`rna_best_pos`, `blast_supported_pos`, `current_start_pos`)とGCAの開始位置
を照合するのみに留める。

ただし、ゲノムFASTA(GCF/GCA完全一致を確認済み)だけで判定できる範囲の
「GCA提案の妥当性チェック」は実施する(Step 3-Bで言及されている
「GCAの提案位置がストップコドン間近」の判定に相当):

  - GCA提案位置(上流延長方向)から現行開始位置までの区間に、読み枠内で
    ストップコドンが存在しないか(あれば、その延長は生物学的に無効)
  - GCA提案位置のコドン自体がATG/GTG/TTGか

これらは新規のBLAST/RNA-seq実行を伴わない、ゲノム配列だけで完結する
決定的なチェックである。

出力
----
  - full_run/full_results_v2_with_gca.csv (全4,683行 + gca_start / gca_agrees_with / revised_confidence)
  - full_run/gca_review_34_summary.csv (「要確認」34件の再評価詳細)
  - full_run/gca_flagged_341.csv (「支持」だがGCAと不一致の遺伝子、延長妥当性チェック付き)
  - full_run/gca_flagged_850.csv (「中確信度」だがGCAと不一致の遺伝子)
  - full_run/gca_frameshift_genes.csv (差分が3の倍数でない遺伝子、タンパク質配列レベルの確認付き)
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

BASE = Path(__file__).resolve().parent
FULL_RUN = BASE / "full_run"
GCF_GFF = "/home/nyako/projects/ProteinHunter/data/databases/target/methanosarcina_acetivorans/ncbi_dataset/data/GCF_000007345.1/genomic.gff"
GCA_GFF = BASE / "data" / "GCA_000007345.1_ASM734v1_genomic.gff"

sys.path.insert(0, str(BASE))
from gca_gcf_comparison import parse_protein_coding_genes  # noqa: E402

STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODONS = {"ATG", "GTG", "TTG"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def load_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, list] = {}
    current = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = []
            elif current is not None:
                seqs[current].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def codon_at(seq: str, pos: int, strand: str) -> str:
    """posは1-based翻訳開始位置(遺伝子の`current_start_pos`等と同じ規約)。"""
    if strand == "+":
        return seq[pos - 1: pos + 2]
    return revcomp(seq[pos - 3: pos])


def check_upstream_extension_valid(seq: str, current_start: int, candidate_start: int, strand: str) -> Tuple[bool, str]:
    """candidate_start(GCA提案、上流側)からcurrent_start(現行)までの区間に
    読み枠内ストップコドンが無く、candidate_startのコドンがATG/GTG/TTGであるかを確認する。
    """
    start_codon = codon_at(seq, candidate_start, strand)
    if start_codon not in START_CODONS:
        return False, f"候補位置のコドンが非正準({start_codon})"

    step = 3 if strand == "+" else -3
    for p in range(candidate_start, current_start, step):
        codon = codon_at(seq, p, strand)
        if codon in STOP_CODONS:
            return False, f"候補位置と現行位置の間に読み枠内ストップコドン({codon}, pos={p})"
    return True, "妥当(区間にストップコドンなし、開始コドンも正準)"


def translate(nt_seq: str) -> str:
    return "".join(CODON_TABLE.get(nt_seq[i:i + 3], "X") for i in range(0, len(nt_seq) - 2, 3))


def get_orf_protein(seq: str, start: int, end: int, strand: str) -> str:
    """start/endはGFFの`gene`行の座標(翻訳開始点を含む側〜ストップコドンを含む側)。"""
    if strand == "+":
        sub = seq[start - 1:end]
    else:
        sub = revcomp(seq[start - 1:end])
    return translate(sub)


def load_csv_by_locus(path: Path) -> Dict[str, dict]:
    with open(path) as fh:
        return {row["locus_tag"]: row for row in csv.DictReader(fh)}


def to_int_or_none(v: Optional[str]) -> Optional[int]:
    if v is None or v == "":
        return None
    return int(float(v))


def classify_gca_agreement(current: int, rna: Optional[int], blast: Optional[int], gca: Optional[int]) -> str:
    if gca is None:
        return "GCA比較不可"
    if gca == current:
        return "現行(GCF)と一致"
    agrees_blast = blast is not None and gca == blast
    agrees_rna = rna is not None and gca == rna
    if agrees_blast and agrees_rna:
        return "BLAST・RNA-seq両方と一致"
    if agrees_blast:
        return "BLASTと一致"
    if agrees_rna:
        return "RNA-seqと一致"
    return "いずれとも異なる(第3の候補)"


def main() -> None:
    genome = load_fasta(BASE / "data" / "GCF_000007345.1_ASM734v1_genomic.fna")
    seqid = next(iter(genome))
    seq = genome[seqid]
    print(f"ゲノム配列読み込み: {seqid} ({len(seq)}bp)")

    gca_by_locus = load_csv_by_locus(FULL_RUN / "gca_gcf_comparison.csv")
    review34_by_locus = load_csv_by_locus(FULL_RUN / "review_34_final_merged.csv")

    with open(FULL_RUN / "full_results.csv") as fh:
        rows = list(csv.DictReader(fh))

    gcf_genes = parse_protein_coding_genes(GCF_GFF)
    gca_genes = parse_protein_coding_genes(str(GCA_GFF))

    out_rows = []
    review34_out = []
    flagged_341 = []
    flagged_850 = []
    frameshift_rows = []

    for row in rows:
        locus_tag = row["locus_tag"]
        current = to_int_or_none(row["current_start_pos"])
        rna_best = to_int_or_none(row["rna_best_pos"])
        blast_supp = to_int_or_none(row["blast_supported_pos"])
        verdict = row["verdict"]

        gca_row = gca_by_locus.get(locus_tag)
        gca_start = to_int_or_none(gca_row["gca_start"]) if gca_row and gca_row["gca_start"] != "" else None
        frame_shift_flag = gca_row["frame_shift_flag"] if gca_row else ""
        is_frameshift = frame_shift_flag == "True"

        agrees_with = classify_gca_agreement(current, rna_best, blast_supp, gca_start)
        revised_confidence = verdict  # デフォルトは変更なし

        # --- 34件「要確認」: review_34_final_merged.csv の final_classification を出発点に反映(3-A) ---
        r34 = review34_by_locus.get(locus_tag)
        if r34 is not None:
            final_class = r34["final_classification"]
            note = ""
            if gca_start is not None and not is_frameshift:
                if final_class == "判断保留":
                    if agrees_with == "BLAST・RNA-seq両方と一致":
                        revised_confidence = "高確信度相当(GCA・RNA-seq・BLASTの3系統一致、決着提案)"
                        note = "GCAがRNA-seq/BLAST両方の候補と一致(この34件は元々RNA-seqとBLASTが食い違うため異例だが確認された)"
                    elif agrees_with == "BLASTと一致":
                        revised_confidence = "中確信度相当(GCA旧アノテーションがBLAST支持位置と一致、決着提案)"
                        note = "GCAがBLAST支持位置と一致したため決着を提案"
                    elif agrees_with == "RNA-seqと一致":
                        revised_confidence = "中確信度相当(GCA旧アノテーションがRNA-seq候補と一致、決着提案)"
                        note = "GCAがRNA-seq候補と一致したため決着を提案"
                    else:
                        revised_confidence = final_class
                        note = "GCAはRNA-seq/BLASTいずれとも異なる第3候補、または現行と一致。判断保留を維持"
                elif final_class == "アーティファクト":
                    revised_confidence = final_class
                    if agrees_with in ("BLASTと一致", "RNA-seqと一致", "BLAST・RNA-seq両方と一致"):
                        note = f"要再検討(両論併記): GCAは{agrees_with}。既存のアーティファクト判定理由(反復配列/減衰テール等)と矛盾しないか要確認"
                    else:
                        note = "GCAはBLAST/RNA候補と一致せず、アーティファクト判定を維持"
                else:
                    # 真の代替候補 / 中確信度相当(格上げ) はそのまま維持し、GCA一致状況のみ記録
                    revised_confidence = final_class
                    note = f"GCA一致状況: {agrees_with}(既存分類を維持)"
            else:
                revised_confidence = final_class
                note = "GCA比較不可(GCA側にマッチする遺伝子なし)" if gca_start is None else "フレームシフト疑いのため個別確認が必要(gca_frameshift_genes.csv参照)"

            review34_out.append({
                "locus_tag": locus_tag,
                "old_locus_tag": row["old_locus_tag"],
                "final_classification_before": final_class,
                "gca_start": gca_start if gca_start is not None else "",
                "gca_agrees_with": agrees_with,
                "revised_confidence": revised_confidence,
                "note": note,
            })

        # --- 3427件「支持」のうちGCA不一致(341件相当、3-B) ---
        elif verdict == "現行アノテーションを支持" and gca_start is not None and gca_start != current and not is_frameshift:
            gene = gcf_genes[locus_tag]
            if gca_start > current if gene.strand == "+" else gca_start < current:
                # GCF側より短い(GCAはGCFより下流を提案): 既存ORF内部なので自明に妥当
                valid, reason = True, "GCA提案がGCFより下流(短縮方向)のため区間チェック対象外"
            else:
                # GCAが上流(延長方向)を提案: ストップコドン混入・非正準開始コドンをチェック
                valid, reason = check_upstream_extension_valid(seq, current, gca_start, gene.strand)
            has_blast_hit = bool(row["blast_best_species_hit"])
            flagged_341.append({
                "locus_tag": locus_tag,
                "old_locus_tag": row["old_locus_tag"],
                "current_start_pos": current,
                "gca_start": gca_start,
                "gca_extension_valid": valid,
                "gca_extension_check_reason": reason,
                "blast_best_species_hit": row["blast_best_species_hit"],
                "blast_pident": row["blast_pident"],
                "rna_step_ratio": row["rna_step_ratio"],
                "weak_evidence_exclude": (not valid) or (not has_blast_hit),
            })

        # --- 1109件「中確信度」のうちGCA不一致(850件相当、3-C) ---
        elif verdict.startswith("中確信度") and gca_start is not None and gca_start != current and not is_frameshift:
            # make_verdict()と同じ基準: RNA-seqは「ステップ比が閾値3.0を超えた」場合のみ
            # 証拠として採用される(rna_best_posが現行と異なる値でも、閾値未満なら不採用)。
            rna_suggests_change = row["rna_step_ratio_exceeds_threshold"] == "True" and rna_best is not None and rna_best != current
            blast_diff = blast_supp is not None and blast_supp != current
            single_pos = single_source = None
            if rna_suggests_change and not blast_diff:
                single_pos, single_source = rna_best, "RNA-seq"
            elif blast_diff and not rna_suggests_change:
                single_pos, single_source = blast_supp, "BLAST"

            if single_pos is not None and gca_start == single_pos:
                revised_confidence = f"高確信度相当(格上げ): {single_source}+GCA旧アノテーションの2系統一致"
            else:
                revised_confidence = verdict

            flagged_850.append({
                "locus_tag": locus_tag,
                "old_locus_tag": row["old_locus_tag"],
                "current_start_pos": current,
                "single_evidence_source": single_source or "",
                "single_evidence_pos": single_pos if single_pos is not None else "",
                "gca_start": gca_start,
                "gca_agrees_with": agrees_with,
                "revised_confidence": revised_confidence,
            })

        # --- フレームシフト疑い: タンパク質配列レベルの確認 ---
        if is_frameshift and gca_row is not None:
            gcf_gene = gcf_genes[locus_tag]
            gca_gene = gca_genes[gca_row["old_locus_tag"]]
            gcf_protein = get_orf_protein(seq, gcf_gene.start, gcf_gene.end, gcf_gene.strand)
            gca_protein = get_orf_protein(seq, gca_gene.start, gca_gene.end, gca_gene.strand)
            gcf_internal_stop = "*" in gcf_protein[:-1]
            gca_internal_stop = "*" in gca_protein[:-1]
            frameshift_rows.append({
                "locus_tag": locus_tag,
                "old_locus_tag": row["old_locus_tag"],
                "strand": gcf_gene.strand,
                "gcf_start": gcf_gene.start if gcf_gene.strand == "+" else gcf_gene.end,
                "gcf_end": gcf_gene.end if gcf_gene.strand == "+" else gcf_gene.start,
                "gca_start": gca_gene.start if gca_gene.strand == "+" else gca_gene.end,
                "gca_end": gca_gene.end if gca_gene.strand == "+" else gca_gene.start,
                "stop_codon_also_differs": (gcf_gene.end != gca_gene.end) if gcf_gene.strand == "+" else (gcf_gene.start != gca_gene.start),
                "diff_bp": gca_row["diff_bp"],
                "diff_codons": gca_row["diff_codons"],
                "gcf_protein_len_aa": len(gcf_protein) - 1,
                "gca_protein_len_aa": len(gca_protein) - 1,
                "gcf_protein_first20aa": gcf_protein[:20],
                "gca_protein_first20aa": gca_protein[:20],
                "gcf_has_unexpected_internal_stop": gcf_internal_stop,
                "gca_has_unexpected_internal_stop": gca_internal_stop,
            })

        out_row = dict(row)
        out_row["gca_start"] = gca_start if gca_start is not None else ""
        out_row["gca_agrees_with"] = agrees_with
        out_row["revised_confidence"] = revised_confidence
        out_rows.append(out_row)

    fieldnames = list(rows[0].keys()) + ["gca_start", "gca_agrees_with", "revised_confidence"]
    with open(FULL_RUN / "full_results_v2_with_gca.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(out_rows)

    with open(FULL_RUN / "gca_review_34_summary.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["locus_tag", "old_locus_tag", "final_classification_before", "gca_start", "gca_agrees_with", "revised_confidence", "note"])
        writer.writeheader()
        writer.writerows(review34_out)

    with open(FULL_RUN / "gca_flagged_341.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["locus_tag", "old_locus_tag", "current_start_pos", "gca_start", "gca_extension_valid", "gca_extension_check_reason", "blast_best_species_hit", "blast_pident", "rna_step_ratio", "weak_evidence_exclude"])
        writer.writeheader()
        writer.writerows(flagged_341)

    with open(FULL_RUN / "gca_flagged_850.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["locus_tag", "old_locus_tag", "current_start_pos", "single_evidence_source", "single_evidence_pos", "gca_start", "gca_agrees_with", "revised_confidence"])
        writer.writeheader()
        writer.writerows(flagged_850)

    with open(FULL_RUN / "gca_frameshift_genes.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["locus_tag", "old_locus_tag", "strand", "gcf_start", "gcf_end", "gca_start", "gca_end", "stop_codon_also_differs", "diff_bp", "diff_codons", "gcf_protein_len_aa", "gca_protein_len_aa", "gcf_protein_first20aa", "gca_protein_first20aa", "gcf_has_unexpected_internal_stop", "gca_has_unexpected_internal_stop"])
        writer.writeheader()
        writer.writerows(frameshift_rows)

    upgraded_34 = sum(1 for r in review34_out if r["revised_confidence"] != r["final_classification_before"])
    upgraded_850 = sum(1 for r in flagged_850 if r["revised_confidence"].startswith("高確信度相当"))
    excluded_341 = sum(1 for r in flagged_341 if r["weak_evidence_exclude"])

    print(f"全{len(out_rows)}行 -> full_results_v2_with_gca.csv")
    print(f"34件再評価 -> gca_review_34_summary.csv ({len(review34_out)}件中 {upgraded_34}件で分類変更提案)")
    print(f"341件相当フラグ -> gca_flagged_341.csv ({len(flagged_341)}件中 {excluded_341}件が根拠薄弱で除外対象)")
    print(f"850件相当フラグ -> gca_flagged_850.csv ({len(flagged_850)}件中 {upgraded_850}件で高確信度相当へ格上げ提案)")
    print(f"フレームシフト疑い -> gca_frameshift_genes.csv ({len(frameshift_rows)}件)")


if __name__ == "__main__":
    main()

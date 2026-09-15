#!/usr/bin/env python3
"""
finalize_reports.py

目的
----
開始コドン検証プロジェクトの最終成果物整理(`final_cleanup_instructions.md`)。
`full_run/`配下の中間ファイル群を統合し、`FinalReports/`配下に「最終的に
参照すべきデータ一式」として配置する。統合結果を書き出すのみで、元ファイル
の削除は行わない(削除は本スクリプトの実行後、行数確認を経て別途行う)。

行う統合
--------
1. `full_run/full_results_v2_with_gca.csv` -> `FinalReports/full_results_final.csv`
   (`FinalReports/full_results.csv`を置き換える。GCA統合前のfull_results.csv
   を完全に上位互換していることを事前に確認済み)
2. `FinalReports/review_34_final_merged.csv`(3段階の深掘り結果)+
   `full_run/gca_review_34_summary.csv`(GCA証拠の反映)
   -> `FinalReports/review_34_final.csv`
   (locus_tagで突き合わせ、GCA側の列を末尾に追加列として残す。
   `final_classification_before`はmerged側の`final_classification`と
   重複するが、GCA証拠側の情報として残す)
3. `full_run/gca_flagged_341.csv`(336行)+
   `full_run/gca_flagged_341_reevaluated.csv`(207行、新規計算済み)
   -> `FinalReports/gca_flagged_support_final.csv`(336行)
   (新規計算されなかった129件には`weak_evidence_exclude_reason`列に
   除外理由を記録)
4. `full_run/gca_flagged_850.csv`(849行)+
   `full_run/gca_flagged_850_reevaluated.csv`(119行、新規計算済み)
   -> `FinalReports/gca_flagged_medium_final.csv`(849行)
5. `full_run/gca_frameshift_genes.csv` -> `FinalReports/gca_frameshift_final.csv`(そのままコピー)
6. `full_run/gca_gcf_comparison.csv` -> `FinalReports/gca_gcf_comparison.csv`(そのままコピー)
"""

from __future__ import annotations

import csv
import shutil
from pathlib import Path
from typing import Dict, List

BASE = Path(__file__).resolve().parent
FULL_RUN = BASE / "full_run"
FINAL = BASE / "FinalReports"


def read_rows(path: Path) -> List[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh))


def read_fieldnames(path: Path) -> List[str]:
    with open(path) as fh:
        return next(csv.reader(fh))


def write_rows(path: Path, rows: List[dict], fieldnames: List[str]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def merge_full_results() -> None:
    src = FULL_RUN / "full_results_v2_with_gca.csv"
    dst = FINAL / "full_results_final.csv"
    rows = read_rows(src)
    write_rows(dst, rows, read_fieldnames(src))
    print(f"[1] {dst} ({len(rows)}行)")


def merge_review_34() -> None:
    merged_path = FINAL / "review_34_final_merged.csv"
    gca_path = FULL_RUN / "gca_review_34_summary.csv"
    merged_rows = read_rows(merged_path)
    merged_fields = read_fieldnames(merged_path)
    gca_by_locus: Dict[str, dict] = {r["locus_tag"]: r for r in read_rows(gca_path)}
    gca_extra_fields = ["final_classification_before", "gca_start", "gca_agrees_with", "revised_confidence", "note"]

    out_rows = []
    for row in merged_rows:
        gca_row = gca_by_locus.get(row["locus_tag"], {})
        out = dict(row)
        for f in gca_extra_fields:
            out[f] = gca_row.get(f, "")
        out_rows.append(out)

    fieldnames = merged_fields + gca_extra_fields
    dst = FINAL / "review_34_final.csv"
    write_rows(dst, out_rows, fieldnames)
    print(f"[2] {dst} ({len(out_rows)}行)")


def merge_flagged_341() -> None:
    base_path = FULL_RUN / "gca_flagged_341.csv"
    reeval_path = FULL_RUN / "gca_flagged_341_reevaluated.csv"
    base_rows = read_rows(base_path)
    base_fields = read_fieldnames(base_path)
    reeval_by_locus: Dict[str, dict] = {r["locus_tag"]: r for r in read_rows(reeval_path)}
    reeval_fields = read_fieldnames(reeval_path)
    extra_fields = [f for f in reeval_fields if f not in base_fields]

    out_rows = []
    for row in base_rows:
        out = dict(row)
        reeval_row = reeval_by_locus.get(row["locus_tag"])
        if reeval_row is not None:
            for f in extra_fields:
                out[f] = reeval_row[f]
            out["weak_evidence_exclude_reason"] = ""
        else:
            for f in extra_fields:
                out[f] = ""
            reasons = []
            if row["gca_extension_valid"] == "False":
                reasons.append(row["gca_extension_check_reason"])
            if not row["blast_best_species_hit"]:
                reasons.append("BLASTヒットなし(既存証拠が弱いため新規計算の対象から除外)")
            out["weak_evidence_exclude_reason"] = "; ".join(reasons)
        out_rows.append(out)

    fieldnames = base_fields + extra_fields + ["weak_evidence_exclude_reason"]
    dst = FINAL / "gca_flagged_support_final.csv"
    write_rows(dst, out_rows, fieldnames)
    print(f"[3] {dst} ({len(out_rows)}行, うち新規再評価済み{len(reeval_by_locus)}行)")


def merge_flagged_850() -> None:
    base_path = FULL_RUN / "gca_flagged_850.csv"
    reeval_path = FULL_RUN / "gca_flagged_850_reevaluated.csv"
    base_rows = read_rows(base_path)
    base_fields = read_fieldnames(base_path)
    reeval_by_locus: Dict[str, dict] = {r["locus_tag"]: r for r in read_rows(reeval_path)}
    reeval_fields = read_fieldnames(reeval_path)
    extra_fields = [f for f in reeval_fields if f not in base_fields]

    out_rows = []
    for row in base_rows:
        out = dict(row)
        reeval_row = reeval_by_locus.get(row["locus_tag"])
        for f in extra_fields:
            out[f] = reeval_row[f] if reeval_row is not None else ""
        out_rows.append(out)

    fieldnames = base_fields + extra_fields
    dst = FINAL / "gca_flagged_medium_final.csv"
    write_rows(dst, out_rows, fieldnames)
    print(f"[4] {dst} ({len(out_rows)}行, うち新規再評価済み{len(reeval_by_locus)}行)")


def copy_as_is(src_name: str, dst_name: str, label: str) -> None:
    src = FULL_RUN / src_name
    dst = FINAL / dst_name
    shutil.copyfile(src, dst)
    with open(dst) as fh:
        n = sum(1 for _ in fh) - 1
    print(f"[{label}] {dst} ({n}行)")


def main() -> None:
    FINAL.mkdir(exist_ok=True)
    merge_full_results()
    merge_review_34()
    merge_flagged_341()
    merge_flagged_850()
    copy_as_is("gca_frameshift_genes.csv", "gca_frameshift_final.csv", "5")
    copy_as_is("gca_gcf_comparison.csv", "gca_gcf_comparison.csv", "6")


if __name__ == "__main__":
    main()

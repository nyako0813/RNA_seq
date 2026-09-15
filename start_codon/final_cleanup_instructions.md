# 開始コドン検証プロジェクト — 最終成果物の整理・不要ファイル削除 指示書（Claude Code向け）

## 背景

スモークテスト→フル実行→34件深掘り→GCA/GCF比較まで一連の分析が完了した。プロジェクトの一区切りとして、`start_codon/FinalReports/`を「最終的に参照すべきデータ一式」として整理し、統合済み・重複した中間ファイルは削除する。

## Step 1: 最終データの統合・配置

以下の統合を行い、結果を`start_codon/FinalReports/`に配置する(統合前に各ファイルの行数を確認し、統合後も行数が一致することを確認してから元ファイルを消すこと)。

1. **全遺伝子結果の最終版**: `full_run/full_results_v2_with_gca.csv`(GCA統合済み、4,683行)を`FinalReports/full_results_final.csv`として配置。既存の`FinalReports/full_results.csv`(GCA統合前のもの)は、内容が完全に上位互換されていることを確認した上で置き換える。

2. **「要確認」34件の最終版**: `full_run/review_34_final_merged.csv`(3段階の深掘り結果)と`full_run/gca_review_34_summary.csv`(GCA証拠の反映)を1つのファイルに統合し、`FinalReports/review_34_final.csv`として配置する。列が重複する場合はGCA証拠側の情報を追加列として残すこと。

3. **「支持」336件の最終版**: `full_run/gca_flagged_341.csv`(336行、第1段階の妥当性チェック)と`full_run/gca_flagged_341_reevaluated.csv`(207行、新規BLAST/RNA-seq計算による最終再評価)を1つのファイルに統合し、`FinalReports/gca_flagged_support_final.csv`として配置する。129件の除外分にも「除外理由」列を残し、207件の再評価結果と合わせて336行全てが1ファイルで追える形にすること。

4. **「中確信度」849件の最終版**: 同様に`full_run/gca_flagged_850.csv`と`full_run/gca_flagged_850_reevaluated.csv`を統合し、`FinalReports/gca_flagged_medium_final.csv`として配置する(849行)。

5. **フレームシフト疑い6件**: `full_run/gca_frameshift_genes.csv`(既に両読み枠のBLAST照合結果を含む最終版)をそのまま`FinalReports/gca_frameshift_final.csv`としてコピー。

6. **GCA/GCF比較の生データ**: `full_run/gca_gcf_comparison.csv`(4,115行、開始位置比較の元データ)をそのまま`FinalReports/gca_gcf_comparison.csv`としてコピー(参照用として有用なため残す)。

7. **非正準開始コドン95件**: 既に`FinalReports/noncanonical_support_flagged.csv`として存在済み、変更不要。

8. **Excelレポート**: 既存の`FinalReports/start_codon_full_run_report_v3.xlsx`はGCA統合前の内容(v3)であることに注意。今回は新しいExcelの作成までは求めないが、その旨を報告に明記すること(必要であれば別途依頼する)。

## Step 2: 不要になった中間ファイルの削除

Step 1の統合が完了し、`FinalReports/`側で全データが確認できたら、以下の`start_codon/full_run/`配下の中間ファイルを削除する。

- `full_results.csv`(GCA統合前、`full_results_v2_with_gca.csv`に完全上位互換されている場合)
- `full_results_v2_with_gca.csv`(`FinalReports/full_results_final.csv`に統合済み)
- `review_34_final_merged.csv`
- `gca_review_34_summary.csv`
- `gca_flagged_341.csv` / `gca_flagged_341_reevaluated.csv`
- `gca_flagged_850.csv` / `gca_flagged_850_reevaluated.csv`
- `gca_frameshift_genes.csv`
- `gca_gcf_comparison.csv`
- `review_candidates.csv`(GCA統合前の1,256件データ。`FinalReports/full_results_final.csv`から再現可能なため削除可。ただし、GCA統合後の確信度で「不一致」に該当する行数が1,256件と変わっている可能性があるため、削除前に一言報告すること)

**削除しないもの**:

- `start_codon/mini_test/` 配下一式(スモークテスト用フィクスチャとして温存)
- `start_codon/full_run/orthologs_methanosarcina.faa`とBLAST DBファイル(`.gitignore`対象、再利用のため温存)
- `start_codon/md/`配下のドキュメント一式(`start_codon_analysis_summary.md`, `gca_gcf_comparison_instructions.md`, `gca_gcf_comparison_report.md`等)
- `start_codon/find_true_start_codons.py`, `integrate_gca_evidence.py`, `resume_gca_gcf_step9e.py`等のスクリプト本体とテストファイル
- `mac/`(ゲノム参照データ)

## Step 3: ドキュメントの参照パス更新

`start_codon/md/start_codon_analysis_summary.md`の§7(成果物一覧)に記載されているファイルパスが、Step 1・2の変更後も正しく`FinalReports/`配下を指すように更新する(パスの列挙のみでよく、内容の書き直しは不要)。

## 出力・報告

- 統合前後の行数確認結果(各ファイルで一致していること)
- `FinalReports/`配下の最終的なファイル一覧
- 削除したファイル一覧
- `review_candidates.csv`削除に伴う件数変化の有無について一言
- 通常運用に従いテスト実行・コミット・push

# 開始コドン検証 — 成果物の整理 指示書（Claude Code向け）

## 背景

「要確認」34件について、第1段階(深度・BLAST・隣接遺伝子による分類)→ 第2段階Task A(BLAST弱証拠グループの再評価)→ Task B(遠方カバレッジ再上昇チェック)と3段階の深掘りを行った。この3つの結果を1つの最終テーブルに統合したので、リポジトリ側もこれに合わせて整理する。

## 添付ファイル

- `review_34_final_merged.csv`（この会話に添付。34件全件・全カラムを1行にまとめた最終版）

最終分類の内訳: 真の代替候補6 / 中確信度相当(格上げ)1 / 判断保留22 / アーティファクト5

## Step 1: 統合ファイルの配置

添付の `review_34_final_merged.csv` を `start_codon/full_run/review_34_final_merged.csv` として保存する。

## Step 2: 中間ファイルの削除

以下3つの中間生成物は `review_34_final_merged.csv` に統合済みのため、`full_run/` から削除してよい。

- `review_34_reclassified.csv`（第1段階の出力）
- `review_23_pending_reevaluated.csv`（Task Aの出力）
- `depth_extended_check.csv`（Task Bの出力）

**削除しないもの**（全遺伝子データの一次ソースであり、34件の話とは独立して必要）:

- `full_results.csv`
- `review_candidates.csv`
- `noncanonical_support_flagged.csv`
- `full_run/orthologs_methanosarcina.faa` および対応するBLAST DBファイル
- `mini_test/` 配下一式（スモークテスト用フィクスチャとして温存）

## Step 3: コミット・push

- `git status`で削除・追加対象を確認する。
- コミットメッセージには「34件の3段階分析結果を1ファイルに統合し、中間生成物を整理した」旨を記載する。
- 通常運用に従い、push前に`git fetch && git status`でリモートと同期していることを確認してからpushする。

## 出力・報告

- 削除したファイル一覧とコミットハッシュを報告する。
- `review_34_final_merged.csv`の行数(34)・列数が想定通りであることを確認してから削除を実行する（誤って統合前にデータが欠落していないか一応チェックすること）。

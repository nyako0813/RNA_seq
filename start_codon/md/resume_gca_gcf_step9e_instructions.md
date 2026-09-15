# 開始コドン検証 — GCA/GCF比較 一時停止作業の再開 指示書（Claude Code向け）

## 背景

PC移行(ノートパソコン→デスクトップパソコン)に伴い、`start_codon_analysis_summary.md`の§9-Eで一時停止していたGCA/GCF比較の追加検証を再開する。当時はRNA-seq BAM・オルソログ蛋白質DB・`samtools`がこの分析環境に存在しなかったため未実施だった。

## Step 0: 環境確認(最優先)

以下が新しいデスクトップ環境(`~/projects/RNA-seq`)に揃っているか確認する。

1. RNA-seq BAM: `rnaseq_requantification/bam/`配下に6サンプル(ユーザーによると存在確認済み。ソート・インデックス済みか`samtools quickcheck`等で確認)
2. `samtools`: 以前使っていたmicromamba環境`rnaseq`が残っているか(`which samtools`または`rnaseq_requantification/tools/mamba_root/envs/rnaseq/bin/samtools`のパス確認)。無ければ再構築が必要。
3. オルソログ蛋白質BLAST DB: `start_codon/full_run/orthologs_methanosarcina.faa`とそのBLAST DBファイル(`.gitignore`対象のため、リポジトリには含まれていない。ローカルに残っているか確認。無ければ、以前の手順〈姉妹プロジェクトのMethanosarcina属6種のprotein.faaを連結〉で再構築する)
4. `blastp`/`makeblastdb`が使える状態か

**いずれか欠けている場合は、無理に進めず欠けているものを報告してから指示を仰ぐこと。** 全て揃っていれば以下のStep 1〜3に進む。

## Step 1(最優先): フレームシフト疑い6件のBLAST照合

`start_codon/FinalReports/`または`full_run/gca_frameshift_genes.csv`に記録された6遺伝子(MA_RS02220, MA_RS07005, MA_RS10350, MA_RS11040, MA_RS28585, MA_RS28795)について、GCF側の読み枠で翻訳したタンパク質配列とGCA側の読み枠で翻訳したタンパク質配列、それぞれをオルソログBLAST DBに対して`blastp`検索する。

- どちらの読み枠がオルソログとより高い同一性・より良いアラインメント(qstartが1に近い)で一致するかを比較する。
- 結果を`gca_frameshift_genes.csv`に列を追加する形で更新する(`gcf_frame_best_pident`, `gcf_frame_best_qstart`, `gca_frame_best_pident`, `gca_frame_best_qstart`, `recommended_frame`)。
- 6件それぞれについて、「GCF読み枠が妥当」「GCA読み枠が妥当」「どちらとも判断つかない(両方ヒットが弱い、または拮抗)」のいずれかを結論として付す。

## Step 2: 「支持」336件中、根拠薄弱と判定されなかった207件の新規証拠計算

`gca_flagged_341.csv`の`weak_evidence_exclude=False`の207件について、GCA提案位置(`gca_start`)を新たな候補開始点として:

1. その位置から翻訳したタンパク質配列でBLAST検索し、同一性・qstartを計算(`blast_pident_at_gca`, `blast_qstart_at_gca`)。
2. `samtools depth`でGCA提案位置前後のRNA-seqステップ比を計算(`rna_step_ratio_at_gca`)。
3. 既存パイプラインの判定基準(ステップ比閾値3.0、BLAST同一性等)に照らして、これらの遺伝子を「現行アノテーションを支持のまま」「中確信度相当に格上げ」「高確信度相当に格上げ」のいずれかに再分類する。

出力: `gca_flagged_341_reevaluated.csv`(207件、新規計算列と再分類結果を追加)

## Step 3: 「中確信度」849件中、GCAと不一致の119件の新規証拠計算

`gca_flagged_850.csv`の「GCAが一致せず中確信度を維持」119件について、Step 2と同様の新規計算を行う。GCA提案位置での証拠と、既存の単独証拠(BLASTまたはRNA-seq)を比較し、以下のいずれかに再分類する。

- GCA提案位置が新規証拠でも支持される場合 → 高確信度相当への格上げを検討
- 支持されない場合 → 中確信度のまま

出力: `gca_flagged_850_reevaluated.csv`(119件)

## 出力・報告

- Step 0の環境確認結果(何が揃っていて何が足りなかったか)
- Step 1〜3の結果CSV
- 全体を通じた判定内訳の変化(件数の増減を表で)
- テスト実行・通常のgit運用(コミット・push)を忘れずに
- 完了後、`start_codon_analysis_summary.md`の§9-Eを「完了」に更新する旨をこちらに報告してください(ドキュメント本体の更新はこちらで行います)

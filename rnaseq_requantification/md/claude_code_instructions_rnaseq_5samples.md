# 追加RNA-seqサンプルの発現量定量・統合 指示書

## 背景・目的

`rnaseq_requant_combined.xlsx`（`gene_counts_tpm` シート、既存6サンプル）は、NCBI SRA から取得した *Methanosarcina acetivorans* C2A の RNA-seq リードを再定量したものである。今回、同一BioProject（PRJNA859665 / PRJNA862363, Song et al. 2023, *ISME J*, DOI: 10.1038/s41396-023-01520-y）から、以下5つの追加SRAランをダウンロード済み。これらを**既存サンプル（特にSRR20281611）と完全に同じ手法**で定量し、既存テーブルに統合する。

同一手法で揃える理由：異なるアライナー・アノテーション・TPM計算式を混在させると、read countやTPMの値が既存6サンプルと直接比較できなくなるため。

## 対象ファイル（5サンプル）

| SRA Run | 条件（condition） | biological replicate |
|---|---|---|
| SRR20281612 | acetate_no-respiratory | 2 |
| SRR20281613 | acetate_no-respiratory | 1 |
| SRR20650029 | methanol_no-respiratory | 3 |
| SRR20650030 | methanol_no-respiratory | 2 |
| SRR20650031 | methanol_no-respiratory | 1 |

参考：`SRR20281611`（acetate_no-respiratory, replicate 3）は既存テーブルに定量済み。今回の`SRR20281612`, `SRR20281613`と合わせて、acetate_no-respiratoryはこれで3反復揃う。methanol_humus-respiratoryの3反復目は**不要**（今回のスコープ外）。

## ステップ0：既存パイプラインの特定（最優先・必須）

既存の `rnaseq_requant_combined.xlsx` を作った際に使用したスクリプト・コマンドを、以下を手がかりに特定すること：

- `/home/nyako/python_scripts/` 以下のスクリプト（ファイル名にRNA-seq, requant, TPM, quantといった語を含むものを探索）
- シェルの履歴（`.bash_history` など）に残る `hisat2` / `STAR` / `salmon` / `featureCounts` / `htseq-count` などのコマンド
- 既存テーブルの列構成（`locus_tag`, `old_locus_tag`, `gene_biotype`, `product`, `chr`, `start`, `end`, `strand`, `length_bp` + サンプルごとの `_read_count` / `_TPM`）から逆算できる情報（例：`locus_tag`がRefSeqの `MA_RS xxxxx` 形式 → NCBI RefSeq GCF_000006325.1 のGFFを使用している可能性が高い）

**見つかった場合**：そのスクリプト・パラメータ（トリミング条件、アライナーとオプション、参照ゲノムのバージョン、ストランド指定、TPM計算式）をそのまま流用して以降のステップを実行する。

**見つからない場合**：下記「デフォルトパイプライン」を使用し、実行後の報告で「既存パイプラインは特定できず、デフォルト条件で実行した」旨を必ず明記する（後でユーザーが既存値と突き合わせて検証できるように）。

## デフォルトパイプライン（既存手法が不明な場合のみ）

1. **入力確認**：5サンプルが `_1.fastq.gz` / `_2.fastq.gz`（paired-end）の形で展開済みか確認。`.sra` のままなら `fasterq-dump --split-files <SRR>` で変換する。
2. **QC/トリミング**：`fastp`（デフォルト設定、アダプタ自動検出）
3. **参照ゲノム・アノテーション**：NCBI RefSeq `GCF_000006325.1`（*M. acetivorans* C2A, 染色体 `NC_003552.1`）のゲノムFASTAとGFFを取得。既存テーブルの `chr` 列が `NC_003552.1` であることと一致させる。
4. **マッピング**：`HISAT2`（ゲノムインデックス作成→マッピング）。ストランド指定が不明な場合は unstranded で実行。
5. **定量**：`featureCounts`（`-p --countReadPairs -t gene -g locus_tag`）でCDS/geneごとのread countを算出。
6. **TPM計算**：標準式で算出する。

   ```
   RPK_i = read_count_i / (length_bp_i / 1000)
   TPM_i = RPK_i / sum(RPK) * 1e6
   ```

## ステップ2：既存テーブルへの統合

`rnaseq_requant_combined.xlsx` の `gene_counts_tpm` シートに、既存の列構成を保ったまま以下10列を追加する：

```
SRR20281612_read_count, SRR20281612_TPM,
SRR20281613_read_count, SRR20281613_TPM,
SRR20650029_read_count, SRR20650029_TPM,
SRR20650030_read_count, SRR20650030_TPM,
SRR20650031_read_count, SRR20650031_TPM
```

`locus_tag` をキーに既存4770遺伝子と1対1で対応することを確認し、行の追加・欠落・順序のずれがないか検証する。

## ステップ3：sample_conditionsシートの更新

既存の `sample_conditions` シートに以下5行を追加する：

```csv
SRR20281612,acetate_no-respiratory
SRR20281613,acetate_no-respiratory
SRR20650029,methanol_no-respiratory
SRR20650030,methanol_no-respiratory
SRR20650031,methanol_no-respiratory
```

## ステップ4：検証（必須）

- 各サンプルの総リード数・マッピング率をログに残し、既存6サンプル（特にSRR20281611）と桁が近いか確認する。
- 各サンプルのTPM列の合計が概ね100万に近いか確認する（TPM定義上の妥当性チェック）。
- 同一条件内での再現性チェック：`SRR20281611` / `SRR20281612` / `SRR20281613`（acetate_no-respiratory 3反復）同士、および `SRR20650029` / `030` / `031`（methanol_no-respiratory 3反復）同士でTPMのサンプル間相関係数（Pearson, log2(TPM+1)推奨）を算出し、大きく乖離するサンプルがないか確認する。

## 成果物

- 更新版 `rnaseq_requant_combined.xlsx`（既存構造を維持し、11サンプル分に拡張）
- 使用したコマンド・パラメータ一式（既存パイプラインを流用した場合はその出典、デフォルトを使った場合はその旨を含む）
- 各サンプルのQC/マッピングログ

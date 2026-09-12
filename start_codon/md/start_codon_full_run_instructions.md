# 開始コドン検証パイプライン — フルゲノムバッチ実行 指示書（Claude Code向け）

## 背景・目的

`find_true_start_codons.py` によるミニ5遺伝子スモークテストは実データで正常に完走済み。パイプライン自体は動作確認済みなので、今回は**全ゲノム（アノテーション上の全タンパク質コード遺伝子、偽遺伝子を除く）を対象にフルバッチを実行**し、結果を出力する。

スモークテストで判明した2件の知見を踏まえ、フル実行では以下を必ず反映すること。

1. **偽遺伝子フィルタは適用済み**（`parse_gff_genes()` を `pseudo=true` のGFFレコードを除外するよう修正済み）。フル実行でもこのフィルタが有効であることを確認してから走らせること。
2. **アノテーション上の開始コドン自体が非正準（ATG/GTG/TTG以外）な遺伝子が一定数存在する**（スモークテスト時点で 125/4702、約2.7%）。これらは「現行アノテーションを支持」という判定が出てもその信頼性が本質的に低い。フル実行の出力には、この「現在の開始コドンが非正準かどうか」を独立した列として必ず含めること（例: `current_codon_noncanonical` の真偽値）。

## 入力データ（パスは全てWSL/Ubuntu-24.04環境内）

- **リポジトリ/スクリプト**: `~/projects/RNA-seq/start_codon/find_true_start_codons.py`
- **GFFアノテーション**: `~/projects/RNA-seq/start_codon/data/methanosarcina_acetivorans/ncbi_dataset/data/GCF_000007345.1/genomic.gff`
- **ゲノムFASTA**: 上記と同じ `GCF_000007345.1/` フォルダ内にあるはず（ファイル名は `ls` で確認すること。スモークテストでは `NC_003552.1.fna` 相当のものを使用した）
- **RNA-seq BAM**: `~/projects/RNA-seq/rnaseq_requantification/bam/`
  - 現時点でこのフォルダには **6ファイルのみ存在**（設計時点で想定していた11サンプル中6）。今回のフル実行はこの6ファイルでそのまま進めてよい。残り5サンプルの要否は別途保留事項として最後にまとめて報告し、ブロッカーにはしない。
- **オルソログ蛋白質FASTA（BLAST用）**: **現在存在しない。今回新規に作成する必要がある。**

## Step 1: 環境確認

- `~/projects/RNA-seq/rnaseq_requantification/` 配下（または近傍）に既存の micromamba 環境 `rnaseq` があるはずなので、そこから `samtools` が使えることを確認する（スモークテスト時点のパス例: `rnaseq_requantification/tools/mamba_root/envs/rnaseq/bin/samtools`。移動している可能性があるので `find`/`which` で再確認）。
- `blastp` / `makeblastdb` が使える状態か確認する（同環境内、またはシステム側）。

## Step 2: オルソログFASTAの再構築

スモークテスト時は、姉妹プロジェクト（ProteinInteractionHunter/GenomeHitFinder系）に既にダウンロード済みだった以下6種のMethanosarcina属ゲノム（*M. acetivorans* 自身は除く）の `protein.faa` を連結して作成した。

- *M. mazei*
- *M. barkeri*
- *M. thermophila*
- *M. vacuolata*
- *M. siciliae*
- *M. horonobensis*

今回は当該ファイルが見当たらないため、まず以下を行うこと。

1. ファイルシステム上で上記6種の `protein.faa`（NCBI datasetsフォーマット、`.../ncbi_dataset/data/GCF_*/protein.faa` のような構造のはず）を検索し、所在を特定する。`~/projects/` 配下および、Windows側からマウントされている可能性のある近傍プロジェクト（ProteinInteractionHunter、GenomeHitFinder等）も探索対象に含める。
2. 6種すべて見つかった場合は、*M. acetivorans* を含めずに連結し、`start_codon/full_run/orthologs_methanosarcina.faa` として保存する（スモークテストでは21,734配列だった。近い数になるか確認）。
3. 一部が見つからない場合は、NCBI datasetsから該当種のタンパク質FASTAを新規ダウンロードして補う（ネットワークアクセスが必要）。
4. 最終的な配列数と内訳（種ごとの配列数）を報告すること。

## Step 3: BLAST DB作成

- `start_codon/full_run/` ディレクトリを新設し、その中で `orthologs_methanosarcina.faa` に対して `makeblastdb` を実行する（スモークテストの `mini_test/` と同じ流儀で、FASTAとBLAST DBファイルを同じディレクトリに置き、再利用可能にする）。
- **`mini_test/` 配下は一切変更しない**（今後のスモークテスト用フィクスチャとして温存する）。フル実行用の成果物はすべて新設の `full_run/` 配下に閉じ込め、プロジェクト外や `mini_test/` 外に一時ファイルを残さないこと。

## Step 4: フルバッチ実行

- `find_true_start_codons.py` を、GFF全体（偽遺伝子除外後の全タンパク質コード遺伝子、想定約4700〜4770件）に対して実行する。実際に処理された遺伝子数を必ず報告すること（スモークテストの知見では偽遺伝子除外後の実タンパク質コード遺伝子数は4702件だった。設計時の想定4770件との差異があれば、その理由（偽遺伝子数163件との整合性など）を確認・報告する）。
- 入力: 上記GFF、ゲノムFASTA、6本のBAM、`full_run/orthologs_methanosarcina.faa` のBLAST DB。
- ロジック・閾値はスモークテストで検証済みのものをそのまま使う（アップストリーム候補コドン探索、RNA-seqステップ比の閾値3.0、BLASTのqstart/identity評価、判定カテゴリ: 現行アノテーションを支持 / 中確信度 / 低確信度 など）。閾値やロジックを変更する必要は無い。

## Step 5: 出力

`start_codon/full_run/full_results.csv` に、遺伝子1件につき1行で以下を出力すること。

- locus_tag（新旧両方、GFFに両方あれば）
- 現行アノテーション上の開始コドン位置・配列
- **`current_codon_noncanonical`**（現行開始コドンがATG/GTG/TTG以外なら true）
- 判定（現行アノテーションを支持／中確信度／低確信度 等、既存カテゴリに準拠）
- BLAST qstart・identity（上位ヒット）
- RNA-seqステップ比とその値が閾値3.0を超えているか
- 採用されたアップストリーム候補コドン（あれば、その位置と塩基配列）

## Step 6: 検証

- **回帰確認**: フル実行結果のうち、ミニ5テスト対象だった5遺伝子（MA_4115, MA_2894, MA_1974, MA_0826, MA_0074/radB）を抽出し、スモークテスト時点の結果（`mini_test/mini5_result.csv`）と一致することを確認する。不一致があれば原因を特定してから報告すること。
- 偽遺伝子（`pseudo=true`）が出力に一切含まれていないことを確認する。
- `current_codon_noncanonical` フラグが立っている遺伝子数が、スモークテストの推定（約125件、全体の2.7%程度）と大きく乖離していないか確認する。

## 最終報告に含めるべき内容

- 処理遺伝子総数、判定カテゴリごとの件数
- `current_codon_noncanonical` の件数と、その中で判定が「支持」寄りになっている件数（要注意ケースとして別途リストアップ）
- 中確信度・低確信度と判定された遺伝子の一覧（人手レビュー候補）
- オルソログFASTAの最終構成（種ごとの配列数、どこから取得したか）
- 未解決事項として: RNA-seq BAMが6/11サンプルのみである点（意図的な絞り込みか、残り5サンプルのダウンロード漏れかは未確認のまま）

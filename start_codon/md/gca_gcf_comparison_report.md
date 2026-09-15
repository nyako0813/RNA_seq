# 開始コドン検証 — GCA/GCF比較 結果報告

[指示書](gca_gcf_comparison_instructions.md)に基づき、GCF_000007345.1(現行RefSeq)とGCA_000007345.1(2002年原論文、GenBank)の開始コドン注釈を比較し、既存のRNA-seq・BLASTパイプライン結果(`full_results.csv`)と統合した。

## 実行環境に関する制約(先に明記)

当時のフル実行で使用したRNA-seq BAM・オルソログ蛋白質DB・`samtools`は本環境に現存しない(`blastp`/`makeblastdb`のみ利用可能)。そのため、指示書Step 3-B/3-Cが求める「GCA提案位置での新規BLAST/RNA-seq計算」は**実施していない**(ユーザー確認済み)。代わりに、既存`full_results.csv`の列(`rna_best_pos`・`blast_supported_pos`・`current_start_pos`)とGCAの開始位置を照合するのみに留めた。ただし、ゲノム配列(GCF/GCA完全一致)だけで判定できる範囲——GCA提案位置が正準開始コドン(ATG/GTG/TTG)か、現行位置との間に読み枠内ストップコドンが無いか——は決定的にチェック済み。

## Step 1: ゲノム配列の同一性確認

- `AE010299.1`(GCA)をNCBI datasetsから新規取得し(`start_codon/data/GCA_000007345.1_ASM734v1_genomic.fna`)、既存の`NC_003552.1`(GCF、`start_codon/data/GCF_000007345.1_ASM734v1_genomic.fna`)と1塩基ずつ比較した。
- **結果: 完全一致**(5,751,492bp、差異ゼロ)。以降の座標比較はそのまま信頼できる。

## Step 2: GCA/GCF比較テーブル

スクリプト: [gca_gcf_comparison.py](../gca_gcf_comparison.py) → [full_run/gca_gcf_comparison.csv](../full_run/gca_gcf_comparison.csv)

`gene_biotype=protein_coding`の行のみを対象に、GCFの`old_locus_tag`とGCAの`locus_tag`(共にMA_XXXX形式)でマッチングした。

| 項目 | 件数 |
|---|---|
| GCF protein_coding遺伝子 | 4,702 |
| GCA protein_coding遺伝子 | 4,540 |
| old_locus_tag経由でマッチ | **4,115** |
| うち開始位置が一致 | 2,808 |
| うち開始位置が不一致 | **1,307**(31.8%) |
| うちフレームシフト疑い(差分が3の倍数でない) | **6** |

予備調査の概算(4,125マッチ・1,313不一致・10フレームシフト)とはやや異なる(本結果はより厳密にprotein_coding同士でのみマッチングし、`gene_biotype`が食い違うケース(下記参照)を除外しているため)。

**副次的な発見**: マッチしなかった39件のGCF protein_coding遺伝子は、GCA(2002年)側では`pseudogene`(偽遺伝子)として注釈されていた(例: MA_RS00795→MA_0142)。つまり、2002年から現在までの間に「偽遺伝子→機能遺伝子」への再分類が39件発生していることになる。これは開始コドンの位置とは別種の知見だが、記録に値する。

## Step 3: 既存パイプライン結果との統合・再評価

スクリプト: [integrate_gca_evidence.py](../integrate_gca_evidence.py) → [full_run/full_results_v2_with_gca.csv](../full_run/full_results_v2_with_gca.csv)(全4,683行、`gca_start`・`gca_agrees_with`・`revised_confidence`列を追加)

### GCA一致状況の全体像(4,683遺伝子)

| gca_agrees_with | 件数 |
|---|---|
| 現行(GCF)と一致 | 2,808 |
| GCA比較不可(マッチなし・フレームシフト等) | 568 |
| BLASTと一致 | 441 |
| いずれとも異なる(第3の候補) | 439 |
| BLAST・RNA-seq両方と一致 | 394 |
| RNA-seqと一致 | 33 |

### 3-A. 「要確認」34件の決着

[full_run/gca_review_34_summary.csv](../full_run/gca_review_34_summary.csv)。`review_34_final_merged.csv`の`final_classification`を出発点に、GCAとの一致状況を反映した。

| 最終分類(before) | 内訳 | 提案後の扱い |
|---|---|---|
| 判断保留(22件) | 13件がGCA≡BLAST支持位置 | → 「中確信度相当(決着提案)」に格上げ |
| | 2件がGCA≡RNA-seq候補 | → 「中確信度相当(決着提案)」に格上げ |
| | 7件はGCA比較不可・第3候補・現行一致 | → 判断保留を維持 |
| 真の代替候補(6件) | 全件維持(既にRNA+BLAST一致の高確信度相当) | 変更なし(GCA一致状況を注記) |
| アーティファクト(5件) | **MA_RS07375, MA_RS12405の2件でGCAがBLAST支持位置と一致** | 判定は維持しつつ「要再検討(両論併記)」を明記(予備調査の指摘と完全に一致) |
| 中確信度相当(格上げ,1件) | GCAは第3候補 | 変更なし |

→ **34件中15件で分類変更(格上げ)を提案**、22件の「判断保留」は7件に縮小。

### 3-B. 「現行アノテーションを支持」3,427件のうちGCA不一致

[full_run/gca_flagged_341.csv](../full_run/gca_flagged_341.csv) — 該当**336件**(予備調査の341件と近似)。

各遺伝子についてGCA提案位置の妥当性をゲノム配列のみでチェック(新規BLAST/RNA-seq計算なし):

| 判定 | 件数 |
|---|---|
| 延長として妥当(ストップコドンなし・正準開始コドン) | 282 |
| 無効(区間にストップコドン混入、または非正準開始コドン) | 54 |
| BLASTヒットが元々存在しない(`blast_best_species_hit`空欄) | 77 |
| **除外基準(妥当性なし OR BLASTヒットなし)に該当** | **129** |
| 除外基準に該当せず、人手レビュー対象として残存 | **207** |

除外基準: GCA提案位置が現行位置との間で読み枠内ストップコドンに遮られる(生物学的に無効な延長)、またはそもそもBLASTで近縁種にヒットしていない(独立支持が期待できない)場合を「根拠薄弱」として除外対象とした。残る207件は、RNA-seq/BLASTの閾値をわずかに下回っていただけの可能性があり、新規BAM/オルソログDBが復元でき次第、優先的に再計算する価値がある。

### 3-C. 「中確信度」1,109件のうちGCA不一致

[full_run/gca_flagged_850.csv](../full_run/gca_flagged_850.csv) — 該当**849件**(予備調査の850件とほぼ一致)。

| 項目 | 件数 |
|---|---|
| 単独証拠がBLASTのみ | 845 |
| 単独証拠がRNA-seqのみ | 4 |
| GCAがその単独証拠と一致 → 高確信度相当へ格上げ提案 | **730**(BLAST側728、RNA-seq側2) |
| GCAが一致せず、中確信度を維持 | 119 |

→ **1,109件中730件(全体の66%)で「2系統一致」となり、高確信度相当への格上げを提案できる。**

### フレームシフト疑い6件

[full_run/gca_frameshift_genes.csv](../full_run/gca_frameshift_genes.csv)。予備調査では「10件」とされていたが、`gene_biotype=protein_coding`同士の厳密なマッチングでは**6件**。

タンパク質配列レベルで確認したところ、6件全てで**開始位置だけでなく終止コドン位置(gene終端)もGCA/GCFで異なっており**、翻訳したタンパク質配列も冒頭から完全に別物だった(下表)。ゲノム配列自体は完全一致するため、これは配列の相違ではなく、**同じ`old_locus_tag`に対して2つの注釈が異なるORF/読み枠を割り当てている**ことを意味する(単純な開始コドンの選び直しではない)。

| locus_tag | old_locus_tag | strand | GCF蛋白長(aa) | GCA蛋白長(aa) | GCF冒頭20aa | GCA冒頭20aa |
|---|---|---|---|---|---|---|
| MA_RS02220 | MA_0424 | + | 76 | 87 | MAFNDRGKPFRGKYDNRGDF | MTEENPSGESTTIVEISNLQ |
| MA_RS07005 | MA_1350 | + | 73 | 88 | MGFNDRGNSYRGRDSGRGGR | VKLKWVLMTEETPTGGETVA |
| MA_RS10350 | MA_1988 | - | 80 | 114 | MKILSMARVKLATKLATAKA | MRTKKFDRSKKTLAILLLLC |
| MA_RS11040 | MA_2125 | - | 73 | 72 | MVDASFSLLQTNPQFATLIL | MYSREENGRRFFLIASNKST |
| MA_RS28585 | MA_2101 | + | 55 | 75 | MSEKVLSYPRKGFLNSHGIE | MTAHVRKSVILSSKRVFKFA |
| MA_RS28795 | MA_2789 | - | 74 | 98 | LETFLRIQTCRILAPGKVLV | MEVKLPWRLSSGSRPAGSWR |

いずれも自身の読み枠内で内部ストップコドンは検出されず(両アノテーションとも自己完結したORFとして成立)、両者は互いに矛盾する別々の遺伝子構造を提案していると解釈するのが妥当。個別のBLAST照合(どちらの読み枠がオルソログと一致するか)は新規計算が必要なため、優先確認リストとして別途扱うことを推奨する。なお`MA_RS28795`はGCF側の開始コドンが`ATG`ではなく非正準("L"始まり)であり、この遺伝子は元々`current_codon_noncanonical`の可能性が高い。

## 未解決事項・今後の課題

1. **新規BLAST/RNA-seq計算が未実施**: Step 3-Bの207件・Step 3-Cの119件(GCAと不一致だが単独証拠と不一致)については、GCA提案位置でのBLAST同一性・RNA-seqカバレッジを新規計算すれば、さらに判定を精緻化できる可能性が高い。BAMファイル・オルソログDB・`samtools`の復元が前提。
2. **予備調査との数値差**: マッチ件数(4,115 vs 4,125)・不一致件数(1,307 vs 1,313)・フレームシフト件数(6 vs 10)がいずれも予備調査よりやや少ない。本結果は`gene_biotype=protein_coding`同士の厳密マッチングによるもので、GCA側で偽遺伝子扱いだった39件を意図的に除外している。予備調査がどの基準でマッチングしたかは不明だが、本結果の方が再現可能かつ説明可能。
3. **39件のpseudogene→protein_coding再分類**: 前述の通り、2002年時点で偽遺伝子とされていた遺伝子がRefSeqで機能遺伝子に再分類されている。開始コドンとは別軸の知見だが、これらの遺伝子は現行パイプラインの「支持」判定であっても、そもそも遺伝子性自体が議論の余地がある可能性がある。
4. **フレームシフト疑い6件**: タンパク質配列そのものが別物であり、どちらの読み枠が正しいかはBLAST照合(新規計算)なしには判断できない。

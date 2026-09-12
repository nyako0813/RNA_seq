# 真の開始コドン推定パイプライン 設計メモ

## 背景・目的

NCBIのアノテーション(GFF)はATGを基準に開始コドンを付けているが、
1. オルソログと比較すると開始位置がずれていることがよくある
2. *Methanosarcina acetivorans* を含む古細菌ではATG以外(GTG, TTGなど)も開始コドンになりうる

という2点から、全4770遺伝子について「現行のATGアノテーションは正しいか、
実際にはもっと上流の開始コドンから翻訳が始まっているのではないか」を、
既存のRNA-seqデータ(BAM)と近縁種オルソログとのBLASTP比較を組み合わせて
判定する。

## 使えるデータで何が言えて、何が言えないか(重要な前提)

今回のRNA-seqは断片化・ランダムプライミングを経た通常のshotgunライブラリで、
5'末端特異的な手法(Cappable-seq, differential RNA-seq, Ribo-seqなど)では
ない。そのため、

- **できること**: 「この領域は転写されているか、遺伝子間領域相当か」という
  数十bp単位の判定。上流に伸ばした候補領域が実際にmRNAとして存在するかどうかの
  判定には十分使える。
- **できないこと**: 数塩基〜数十塩基しか離れていない候補コドン同士を、
  RNA-seqカバレッジだけで1塩基単位まで確実に区別すること。

このため、RNA-seqカバレッジによる判定と、オルソログBLASTPによる判定
(保存されたN末端がゲノム上のどこに対応するか)を**独立した2本の証拠**として
組み合わせ、両者が一致した場合のみ「高確信度」とする設計にした。どちらか
一方しか変更を支持しない場合は「中確信度」、両者が異なる候補を指す場合は
「要確認」として、最終判断は人間(Nyakoさん)が行う前提。

## アルゴリズム全体像

各遺伝子について:

1. **候補開始コドンの列挙**
   現行の開始コドンから読み枠を保ったまま上流に走査し、ATG/GTG/TTGが
   出現するたびに候補として記録する。走査は次のいずれかで停止する。
   - 同じ読み枠でストップコドンに当たった(そこより上流は別のORF)
   - 隣接遺伝子(同じ鎖上で直前にある遺伝子)の領域に踏み込む
   - 上限距離(既定300bp)に達した

2. **延長タンパク質配列の作成**
   最も上流の候補から現行の終止コドンまでを翻訳した「延長タンパク質」を
   作る(最初のコドンはATG/GTG/TTGいずれでもMとして翻訳)。これをBLASTPの
   クエリにする。

3. **RNA-seqカバレッジ評価**
   `samtools depth` で各候補コドンの前後(既定30bp窓)の深度を取得し、
   「直後の平均深度 / 直前の平均深度」を立ち上がりスコアとする。全サンプル
   のBAMをプールして深度を合算することで、条件によらない「その領域が
   転写されているかどうか」の頑健な判定にしている(特定条件だけで発現する
   遺伝子を見落とさないため)。

4. **オルソログBLASTP評価**
   延長タンパク質を近縁種のタンパク質データベースに対してBLASTPし、
   最良ヒットのクエリ側アラインメント開始位置(qstart)が、延長タンパク質の
   中のどの候補に最も近いかを求める。qstartが延長配列の先頭(=最も上流の候補)
   に近ければ、上流への延長を支持する強い証拠になる。

5. **統合判定**
   - RNA-seqとBLASTPが同じ上流候補を支持 → **高確信度**で上流延長を提案
   - どちらか一方のみが上流候補を支持 → **中確信度**(要目視確認)
   - 両者が異なる候補を支持 → **要確認**(自動判定は保留)
   - どちらも現行アノテーションを支持 → **現行アノテーションを支持**

## 必要な準備

### 1. ソフトウェア
- `samtools`(既存パイプラインですでに使用済みのはず)
- NCBI BLAST+ (`makeblastdb`, `blastp`)
  - WSL上で未導入なら: `sudo apt install ncbi-blast+` または
    `conda install -c bioconda blast`

### 2. オルソログ用プロテインFASTA
近縁種(推奨: *Methanosarcina mazei* Go1, *Methanosarcina barkeri* Fusaro,
*Methanosarcina thermophila* など)のRefSeqタンパク質FASTAを取得し、
1つのファイルに連結する。NCBI datasetsツールを使う場合の例:

```bash
datasets download genome accession GCF_000006605.1 GCF_000969205.1 \
    --include protein
# 展開して protein.faa をすべて cat で連結し orthologs.faa を作る
cat */protein.faa > orthologs.faa
```

自分が使い慣れた方法(NCBI Webから手動ダウンロードでも可)で問題ない。
複数種混ぜるほどヒット率が上がるが、種数はいくつでもよい。

### 3. BAMファイル
既存パイプラインで作成済みの各サンプルのソート済み・インデックス済みBAM
(`samtools index` 済み)をそのまま使う。全11サンプルをまとめて渡せば
スクリプト内でプールして深度を合算する。

## 実行例

```bash
python find_true_start_codons.py \
    --genome genome/NC_003552.1.fna \
    --gff genome/GCF_000007345.1/genomic.gff \
    --bam bam/SRR20281608.sorted.bam bam/SRR20281609.sorted.bam \
          bam/SRR20281610.sorted.bam bam/SRR20281611.sorted.bam \
          bam/SRR20280855.sorted.bam bam/SRR20280860.sorted.bam \
          bam/SRR20281612.sorted.bam bam/SRR20281613.sorted.bam \
          bam/SRR20650029.sorted.bam bam/SRR20650030.sorted.bam \
          bam/SRR20650031.sorted.bam \
    --orthologs orthologs.faa \
    --out start_codon_reannotation.csv
```

4770遺伝子分のCDS抽出・カバレッジ取得・BLASTPは、ローカルBLASTであれば
数分〜十数分程度で終わる見込み(BAMアクセスがボトルネックになりやすいので、
BAMはSSD上に置くことを推奨)。

## 出力(CSV)の列

| 列名 | 内容 |
|---|---|
| `locus_tag` | 遺伝子ID |
| `current_start_pos` / `current_codon` | 現行アノテーションの開始位置・コドン |
| `n_candidates_upstream` | 上流に見つかった候補コドン数 |
| `upstream_limit_reason` | 探索がどこで止まったか(`stop_codon`/`neighbor_gene`/`max_distance`/`contig_end`) |
| `rna_best_pos` / `rna_best_codon` / `rna_step_ratio` | RNA-seqが支持する候補とその立ち上がりスコア |
| `blast_best_species_hit` / `blast_qstart_aa` / `blast_pident` | BLASTPの最良ヒットとその情報 |
| `blast_supported_pos` / `blast_supported_codon` | BLASTPが支持する候補 |
| `proposed_start_pos` / `proposed_codon` | 統合判定による提案開始位置(RNA-seqとBLASTが一致する場合はその位置) |
| `verdict` | 最終判定(「現行アノテーションを支持」「高確信度」「中確信度」「要確認」) |

## 使い方の推奨フロー

1. まず `verdict` が「高確信度」の行だけ抽出して確認する。ここが最も
   確実に見直す価値のある遺伝子。
2. 次に「中確信度」の行を目視で確認する(IGV等でBAMとGFFを並べて見る
   のが確実)。
3. MA_4115パスウェイの候補遺伝子(MA_4115, MA_2894, MA_1974, MA_0826,
   MA_0074)は件数が少ないので、`verdict`に関わらず個別にIGVで目視確認
   することを推奨する。

## 既知の制約・注意点

- ライブラリがstrand-specific(鎖特異的)かどうかが不明なため、今回の
  カバレッジ取得は鎖を区別せず両鎖の合計depthを使っている。もし対向鎖
  (アンチセンス転写やお隣の遺伝子)の発現が強い遺伝子座では、カバレッジ
  スコアが実際より高く/低く出る可能性がある。ライブラリがstrand-specific
  だと分かっている場合は教えてほしい。`samtools depth`は鎖情報を区別
  しないため、対応する場合は `samtools view -f/-F` でstrand別に分けて
  カウントする改修が必要になる。
- `--max-upstream` (既定300bp)は遺伝子間領域が短い原核生物ゲノムでは
  妥当な値だが、極端に長いIGRを持つ遺伝子座では変更が必要な場合がある。
- BLASTPの`-max_target_seqs 1`は「最良1ヒットのみ」を意味しない場合が
  ある(BLAST+の既知の仕様上の癖)。より厳密にしたい場合は
  `-max_target_seqs 5`程度にして`parse_best_blast_hits`側で
  bitscore最大を選ぶ現在の実装のままで問題なく動作する。

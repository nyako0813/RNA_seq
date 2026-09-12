#!/usr/bin/env python3
"""
find_true_start_codons.py

目的
----
GFF上でATG基準にアノテーションされている各遺伝子について、
  (1) 上流に存在しうる代替開始コドン(ATG/GTG/TTG)を読み枠内で列挙し、
  (2) RNA-seqのカバレッジ(BAM)がどこで「遺伝子間領域相当」から
      「転写産物相当」に立ち上がるかを見て候補を評価し、
  (3) 近縁種オルソログとのBLASTPアラインメントで、保存されたN末端が
      ゲノム上のどの候補に対応するかを見て候補を評価し、
(2)と(3)の証拠を統合して「現行アノテーションを支持する」か
「より上流の開始コドンを提案する」かを遺伝子ごとに判定する。

想定環境
--------
- samtools (`samtools depth`) が PATH に通っていること
- NCBI BLAST+ (`makeblastdb`, `blastp`) が PATH に通っていること
- Python 3.8+ (追加ライブラリ不要。標準ライブラリのみで動作します)

使い方の概要
------------
    python find_true_start_codons.py \\
        --genome genome/NC_003552.1.fna \\
        --gff genome/GCF_000007345.1/genomic.gff \\
        --bam counts/bam/SRR20281611.sorted.bam counts/bam/SRR20281612.sorted.bam ... \\
        --orthologs orthologs/related_species_proteins.faa \\
        --out start_codon_reannotation.csv

`--orthologs` を省略した場合はBLASTのステップをスキップし、RNA-seq
カバレッジの証拠のみで判定します(詳細は設計メモを参照)。

このスクリプトは「初心者でも中身を追いながら手を入れられる」ことを
優先し、あえて1ファイル・素朴な実装にしています。関数ごとに独立して
テストできるようになっているので、まずは candidate 抽出だけ動かす、
次にカバレッジだけ足す、という順番で試すこともできます。
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
import urllib.parse
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# 定数
# ---------------------------------------------------------------------------

# 古細菌(Archaea)で報告されている主な開始コドン。頻度順ではなく、
# 単純にATG以外の候補として何を許容するかのリストです。
ALT_START_CODONS = {"ATG", "GTG", "TTG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}

# 標準コドン表(翻訳表11はStopの構成が標準表と同じなので、開始コドンの
# 特別扱い以外はこれで問題ありません)。
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

# 上流探索の上限(bp)。ここまで遡ってもin-frameのストップコドンが
# 見つからない場合はそこで打ち切る(=隣の遺伝子に食い込まない範囲で)。
MAX_UPSTREAM_SEARCH = 300
# カバレッジのステップ判定に使う窓幅(bp)
COVERAGE_WINDOW = 30


# ---------------------------------------------------------------------------
# データ構造
# ---------------------------------------------------------------------------

@dataclass
class Gene:
    locus_tag: str
    seqid: str
    start: int  # 1-based, GFF準拠 (CDSの最初の塩基、strand=+ ならATGのA)
    end: int    # 1-based, GFF準拠 (CDSの最後の塩基、strand=+ ならストップコドンの最後)
    strand: str  # "+" or "-"
    old_locus_tag: Optional[str] = None  # 例: "MA_4115"(GFFのold_locus_tag属性由来)


@dataclass
class Candidate:
    # ゲノム上の座標(1-based, 開始コドンの最初の塩基)
    genomic_pos: int
    codon: str
    # 現行アノテーションの開始コドンから数えて何コドン上流か(0 = 現行そのもの)
    codons_upstream_of_current: int
    # 「延長タンパク質」(下記参照)の中でのアミノ酸オフセット(1-based)
    aa_offset_in_extended: int


@dataclass
class GeneCandidates:
    gene: Gene
    candidates: List[Candidate]  # 上流側から順、最後が現行の開始コドン
    extended_aa_seq: str  # 一番上流の候補から翻訳したタンパク質配列(BLAST用)
    upstream_limit_reason: str  # "stop_codon" | "neighbor_gene" | "max_distance"


# ---------------------------------------------------------------------------
# FASTA / GFF パース(外部ライブラリ非依存の素朴な実装)
# ---------------------------------------------------------------------------

def load_fasta(path: str) -> Dict[str, str]:
    """マルチFASTAを {配列ID: 配列文字列(大文字)} で返す。"""
    seqs: Dict[str, List[str]] = {}
    current_id = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                seqs[current_id] = []
            else:
                if current_id is None:
                    raise ValueError(f"FASTAの形式が不正です: {path}")
                seqs[current_id].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def parse_old_locus_tag(raw: Optional[str]) -> Optional[str]:
    """GFFのold_locus_tag属性値(例: "MA4115%2CMA_4115")から、
    "MA_XXXX"形式(アンダースコア付き)のタグを取り出す。無ければ元の値をそのまま返す。
    """
    if not raw:
        return None
    decoded = urllib.parse.unquote(raw)
    candidates = [v.strip() for v in decoded.split(",") if v.strip()]
    for c in candidates:
        if re.fullmatch(r"MA_\d+", c):
            return c
    return decoded


def parse_gff_genes(path: str, feature_type: str = "CDS") -> List[Gene]:
    """GFF3から遺伝子(既定ではCDS行)を抽出する。

    NCBI RefSeqのGFF3を想定し、属性欄から locus_tag を取り出す。
    locus_tag が無い場合は ID/Name で代用する。

    old_locus_tag属性は"gene"行にのみ付与され"CDS"行には無いため、
    先に"gene"行だけを走査してlocus_tag毎のold_locus_tagを引けるようにしておく。
    """
    old_locus_tag_by_locus: Dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9 or cols[2] != "gene":
                continue
            attr_dict = dict(kv.split("=", 1) for kv in cols[8].split(";") if "=" in kv)
            locus_tag = attr_dict.get("locus_tag")
            old_locus_tag = parse_old_locus_tag(attr_dict.get("old_locus_tag"))
            if locus_tag and old_locus_tag:
                old_locus_tag_by_locus[locus_tag] = old_locus_tag

    genes: List[Gene] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, frame, attrs = cols
            if ftype != feature_type:
                continue
            attr_dict = {}
            for kv in attrs.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attr_dict[k] = v
            if attr_dict.get("pseudo") == "true":
                # 偽遺伝子(NCBIが不完全/機能喪失と判定した断片)は開始コドン
                # 再推定の対象にならない(N末端が失われているなど、記録された
                # 座標が本来の開始コドンに対応しない)ためスキップする。
                continue
            locus_tag = attr_dict.get("locus_tag") or attr_dict.get("gene") or attr_dict.get("ID")
            if locus_tag is None:
                continue
            genes.append(Gene(
                locus_tag=locus_tag,
                seqid=seqid,
                start=int(start),
                end=int(end),
                strand=strand,
                old_locus_tag=old_locus_tag_by_locus.get(locus_tag),
            ))

    # 同じlocus_tagが複数CDS行を持つ場合(例: プログラム化リボソームフレームシフトで
    # ORF-A/ORF-Bに分かれて注釈されているIS因子トランスポゼース)は、本スクリプトの
    # 「1遺伝子=1連続CDS」という前提が成り立たない。無理に結合すると誤った上流探索・
    # 誤った延長タンパク質になるため、そうしたlocus_tagは丸ごと除外する。
    counts: Dict[str, int] = defaultdict(int)
    for g in genes:
        counts[g.locus_tag] += 1
    duplicated = sorted(lt for lt, n in counts.items() if n > 1)
    if duplicated:
        print(
            f"  警告: {len(duplicated)}件のlocus_tagが複数CDS行を持つため除外します "
            f"(プログラム化フレームシフト等、1連続CDSを前提とする本スクリプトでは扱えません): "
            f"{', '.join(duplicated)}",
            file=sys.stderr,
        )
    return [g for g in genes if counts[g.locus_tag] == 1]


# ---------------------------------------------------------------------------
# 配列ユーティリティ
# ---------------------------------------------------------------------------

def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def get_codon(genome_seq: str, seqid_seq: str, pos_1based: int, strand: str) -> Optional[str]:
    """strand基準で pos_1based を先頭とする3塩基コドンを返す。
    pos_1based は「+鎖ならその座標から3'方向に3塩基」「-鎖ならその座標から
    5'方向(ゲノム座標としては小さくなる方向)に3塩基を取ってから逆相補」を意味する。
    範囲外なら None。
    """
    if strand == "+":
        s = pos_1based - 1
        e = s + 3
        if s < 0 or e > len(seqid_seq):
            return None
        return seqid_seq[s:e]
    else:
        e = pos_1based
        s = e - 3
        if s < 0 or e > len(seqid_seq):
            return None
        return revcomp(seqid_seq[s:e])


def step_genomic_pos(pos_1based: int, strand: str, n_codons: int) -> int:
    """strand基準で n_codons コドン分「上流」に移動した後のゲノム座標を返す
    (n_codonsは正の整数、上流方向への移動)。
    """
    if strand == "+":
        return pos_1based - 3 * n_codons
    else:
        return pos_1based + 3 * n_codons


def translate_cds(seq: str) -> str:
    """3の倍数長のCDS配列(最初のコドンは開始コドンとして扱いM固定)をアミノ酸に翻訳する。"""
    if len(seq) < 3:
        return ""
    aa = ["M"]  # 開始コドンはATG/GTG/TTGいずれでもMとして翻訳するのが慣例
    for i in range(3, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa)


# ---------------------------------------------------------------------------
# ステップ1: 候補開始コドンの列挙
# ---------------------------------------------------------------------------

def find_candidates(
    gene: Gene,
    genome: Dict[str, str],
    all_genes_by_seqid_strand: Dict[Tuple[str, str], List[Gene]],
    max_upstream: int = MAX_UPSTREAM_SEARCH,
) -> GeneCandidates:
    """現行の開始コドンから上流に読み枠を保ったまま遡り、
    ATG/GTG/TTG候補とその上限(ストップコドン/隣接遺伝子/距離上限)を求める。
    """
    seq = genome[gene.seqid]
    current_start_pos = gene.start if gene.strand == "+" else gene.end

    # 隣接遺伝子(同じseqid・同じstrand上で自分より上流側にある最も近い遺伝子)の
    # 終端を求め、そこを越えて探索しないようにする。
    neighbor_limit = None
    same_strand_genes = all_genes_by_seqid_strand.get((gene.seqid, gene.strand), [])
    if gene.strand == "+":
        upstream_ends = [g.end for g in same_strand_genes if g.end < gene.start and g.locus_tag != gene.locus_tag]
        if upstream_ends:
            neighbor_limit = max(upstream_ends)
    else:
        upstream_ends = [g.start for g in same_strand_genes if g.start > gene.end and g.locus_tag != gene.locus_tag]
        if upstream_ends:
            neighbor_limit = min(upstream_ends)

    candidates: List[Candidate] = []
    reason = "max_distance"
    n_codons = 0
    while True:
        n_codons += 1
        pos = step_genomic_pos(current_start_pos, gene.strand, n_codons)

        # 距離上限チェック
        if 3 * n_codons > max_upstream:
            reason = "max_distance"
            break

        # 隣接遺伝子チェック(その遺伝子の領域に踏み込む前に停止)
        if neighbor_limit is not None:
            # 「+」鎖: このコドンの最小座標(pos)が隣接遺伝子の終端以下なら重なる
            # 「-」鎖: このコドンの最大座標(pos)が隣接遺伝子の開始以上なら重なる
            if gene.strand == "+" and pos <= neighbor_limit:
                reason = "neighbor_gene"
                break
            if gene.strand == "-" and pos >= neighbor_limit:
                reason = "neighbor_gene"
                break

        codon = get_codon(genome[gene.seqid], seq, pos, gene.strand)
        if codon is None:
            reason = "contig_end"
            break
        if codon in STOP_CODONS:
            reason = "stop_codon"
            break
        if codon in ALT_START_CODONS:
            candidates.append(Candidate(
                genomic_pos=pos,
                codon=codon,
                codons_upstream_of_current=n_codons,
                aa_offset_in_extended=0,  # 後で埋める
            ))

    # 現行の開始コドンも候補として末尾に追加(codons_upstream_of_current=0)
    current_codon = get_codon(genome[gene.seqid], seq, current_start_pos, gene.strand)
    candidates.append(Candidate(
        genomic_pos=current_start_pos,
        codon=current_codon or "???",
        codons_upstream_of_current=0,
        aa_offset_in_extended=0,
    ))

    # 上流順(codons_upstream_of_current が大きい順)に並べ直す
    candidates.sort(key=lambda c: -c.codons_upstream_of_current)

    # 延長タンパク質配列を作る: 一番上流の候補から現行遺伝子のストップコドンまで
    most_upstream = candidates[0]
    if gene.strand == "+":
        ext_start = most_upstream.genomic_pos
        ext_end = gene.end
        ext_dna = seq[ext_start - 1:ext_end]
    else:
        ext_start = gene.start
        ext_end = most_upstream.genomic_pos
        ext_dna = revcomp(seq[ext_start - 1:ext_end])
    ext_aa = translate_cds(ext_dna)

    # 各候補のアミノ酸オフセットを計算(1-based, 延長配列内での位置)
    for c in candidates:
        c.aa_offset_in_extended = c.codons_upstream_of_current - most_upstream.codons_upstream_of_current + 1

    return GeneCandidates(
        gene=gene,
        candidates=candidates,
        extended_aa_seq=ext_aa,
        upstream_limit_reason=reason,
    )


def build_genes_by_seqid_strand(genes: List[Gene]) -> Dict[Tuple[str, str], List[Gene]]:
    d: Dict[Tuple[str, str], List[Gene]] = defaultdict(list)
    for g in genes:
        d[(g.seqid, g.strand)].append(g)
    return d


# ---------------------------------------------------------------------------
# ステップ2: RNA-seqカバレッジによる評価
# ---------------------------------------------------------------------------

def get_pooled_depth(
    bam_paths: List[str],
    seqid: str,
    region_start: int,
    region_end: int,
    samtools_bin: str = "samtools",
) -> Dict[int, int]:
    """samtools depth を使って複数BAMのdepthを合算し、{ゲノム座標(1-based): 深度} を返す。"""
    region = f"{seqid}:{region_start}-{region_end}"
    cmd = [samtools_bin, "depth", "-a", "-r", region] + bam_paths
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    depth: Dict[int, int] = {}
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        pos = int(parts[1])
        total = sum(int(x) for x in parts[2:])
        depth[pos] = total
    return depth


def mean_depth(depth: Dict[int, int], start: int, end: int) -> float:
    vals = [depth.get(p, 0) for p in range(start, end + 1)]
    if not vals:
        return 0.0
    return sum(vals) / len(vals)


def score_candidates_by_coverage(
    gc: GeneCandidates,
    depth: Dict[int, int],
    window: int = COVERAGE_WINDOW,
) -> Dict[int, float]:
    """各候補(genomic_posをキー)について、
    (直後 window bp の平均depth + 1) / (直前 window bp の平均depth + 1)
    を「立ち上がりの鋭さ」スコアとして返す。値が大きいほど
    「その手前は転写されておらず、そこから発現が始まっている」ことを示唆する。
    """
    scores: Dict[int, float] = {}
    strand = gc.gene.strand
    for c in gc.candidates:
        if strand == "+":
            up_start, up_end = c.genomic_pos - window, c.genomic_pos - 1
            down_start, down_end = c.genomic_pos, c.genomic_pos + window - 1
        else:
            up_start, up_end = c.genomic_pos + 1, c.genomic_pos + window
            down_start, down_end = c.genomic_pos - window + 1, c.genomic_pos
        up_mean = mean_depth(depth, up_start, up_end)
        down_mean = mean_depth(depth, down_start, down_end)
        scores[c.genomic_pos] = (down_mean + 1.0) / (up_mean + 1.0)
    return scores


# ---------------------------------------------------------------------------
# ステップ3: オルソログBLASTPによる評価
# ---------------------------------------------------------------------------

def write_extended_fasta(gene_candidates: List[GeneCandidates], out_path: str) -> None:
    with open(out_path, "w") as fh:
        for gc in gene_candidates:
            if gc.extended_aa_seq:
                fh.write(f">{gc.gene.locus_tag}\n{gc.extended_aa_seq}\n")


def run_blast(
    query_fasta: str,
    ortholog_fasta: str,
    out_tsv: str,
    makeblastdb_bin: str = "makeblastdb",
    blastp_bin: str = "blastp",
    threads: int = 4,
    evalue: float = 1e-5,
) -> None:
    db_marker = ortholog_fasta + ".psq"
    if not Path(db_marker).exists():
        subprocess.run(
            [makeblastdb_bin, "-in", ortholog_fasta, "-dbtype", "prot"],
            check=True,
        )
    subprocess.run(
        [
            blastp_bin,
            "-query", query_fasta,
            "-db", ortholog_fasta,
            "-out", out_tsv,
            "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
            "-evalue", str(evalue),
            "-max_target_seqs", "1",
            "-num_threads", str(threads),
        ],
        check=True,
    )


def parse_best_blast_hits(tsv_path: str) -> Dict[str, dict]:
    """qseqid ごとにbitscore最大の行を採用して返す。"""
    best: Dict[str, dict] = {}
    with open(tsv_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 10:
                continue
            qseqid, sseqid, pident, length, qstart, qend, sstart, send, evalue, bitscore = parts
            rec = {
                "sseqid": sseqid,
                "pident": float(pident),
                "length": int(length),
                "qstart": int(qstart),
                "qend": int(qend),
                "evalue": float(evalue),
                "bitscore": float(bitscore),
            }
            if qseqid not in best or rec["bitscore"] > best[qseqid]["bitscore"]:
                best[qseqid] = rec
    return best


def blast_supported_candidate(gc: GeneCandidates, best_hit: Optional[dict]) -> Optional[Candidate]:
    """BLASTのqstart(延長タンパク質中のアミノ酸位置)に最も近い候補を返す。"""
    if best_hit is None:
        return None
    qstart = best_hit["qstart"]
    return min(gc.candidates, key=lambda c: abs(c.aa_offset_in_extended - qstart))


# ---------------------------------------------------------------------------
# ステップ4: 統合判定
# ---------------------------------------------------------------------------

COVERAGE_STEP_THRESHOLD = 3.0  # このステップ比以上を「立ち上がりあり」とみなす


def make_verdict(
    gc: GeneCandidates,
    coverage_scores: Optional[Dict[int, float]],
    best_blast_hit: Optional[dict],
) -> dict:
    current = next(c for c in gc.candidates if c.codons_upstream_of_current == 0)

    rna_best = None
    rna_best_score = None
    if coverage_scores:
        rna_best_pos = max(coverage_scores, key=lambda p: coverage_scores[p])
        rna_best = next(c for c in gc.candidates if c.genomic_pos == rna_best_pos)
        rna_best_score = coverage_scores[rna_best_pos]

    blast_best = blast_supported_candidate(gc, best_blast_hit) if best_blast_hit else None

    rna_suggests_change = (
        rna_best is not None
        and rna_best.codons_upstream_of_current != 0
        and rna_best_score is not None
        and rna_best_score >= COVERAGE_STEP_THRESHOLD
    )
    blast_suggests_change = blast_best is not None and blast_best.codons_upstream_of_current != 0

    if rna_suggests_change and blast_suggests_change and rna_best.genomic_pos == blast_best.genomic_pos:
        verdict = "高確信度: 上流への延長を提案(RNA-seq・BLAST一致)"
        proposed = rna_best
    elif rna_suggests_change and blast_suggests_change:
        verdict = "要確認: RNA-seqとBLASTが異なる候補を支持"
        proposed = None
    elif rna_suggests_change or blast_suggests_change:
        verdict = "中確信度: 上流への延長の可能性(片方の証拠のみ)"
        proposed = rna_best if rna_suggests_change else blast_best
    else:
        verdict = "現行アノテーションを支持"
        proposed = current

    return {
        "locus_tag": gc.gene.locus_tag,
        "old_locus_tag": gc.gene.old_locus_tag,
        "current_start_pos": current.genomic_pos,
        "current_codon": current.codon,
        "current_codon_noncanonical": current.codon not in ALT_START_CODONS,
        "n_candidates_upstream": len(gc.candidates) - 1,
        "upstream_limit_reason": gc.upstream_limit_reason,
        "rna_best_pos": rna_best.genomic_pos if rna_best else None,
        "rna_best_codon": rna_best.codon if rna_best else None,
        "rna_step_ratio": round(rna_best_score, 2) if rna_best_score is not None else None,
        "rna_step_ratio_exceeds_threshold": (
            rna_best_score is not None and rna_best_score >= COVERAGE_STEP_THRESHOLD
        ),
        "blast_best_species_hit": best_blast_hit["sseqid"] if best_blast_hit else None,
        "blast_qstart_aa": best_blast_hit["qstart"] if best_blast_hit else None,
        "blast_pident": best_blast_hit["pident"] if best_blast_hit else None,
        "blast_supported_pos": blast_best.genomic_pos if blast_best else None,
        "blast_supported_codon": blast_best.codon if blast_best else None,
        "proposed_start_pos": proposed.genomic_pos if proposed else None,
        "proposed_codon": proposed.codon if proposed else None,
        "verdict": verdict,
    }


# ---------------------------------------------------------------------------
# メイン処理
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True, help="ゲノムFASTA")
    ap.add_argument("--gff", required=True, help="GFF3ファイル")
    ap.add_argument("--bam", nargs="+", required=True, help="ソート済みindex済みBAM(複数可、プールして使用)")
    ap.add_argument("--orthologs", help="オルソログのプロテインFASTA(複数種を1ファイルに連結したもの)。省略するとBLASTをスキップ")
    ap.add_argument("--out", default="start_codon_reannotation.csv")
    ap.add_argument("--max-upstream", type=int, default=MAX_UPSTREAM_SEARCH)
    ap.add_argument("--samtools-bin", default="samtools")
    ap.add_argument("--blastp-bin", default="blastp")
    ap.add_argument("--makeblastdb-bin", default="makeblastdb")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--workdir", default=".fsc_work")
    args = ap.parse_args()

    Path(args.workdir).mkdir(exist_ok=True)

    print("[1/5] ゲノム・GFFを読み込み中...", file=sys.stderr)
    genome = load_fasta(args.genome)
    genes = parse_gff_genes(args.gff)
    genes_by_seqid_strand = build_genes_by_seqid_strand(genes)
    print(f"  遺伝子数: {len(genes)}", file=sys.stderr)

    print("[2/5] 候補開始コドンを列挙中...", file=sys.stderr)
    all_gc: List[GeneCandidates] = []
    for g in genes:
        try:
            gc = find_candidates(g, genome, genes_by_seqid_strand, max_upstream=args.max_upstream)
            all_gc.append(gc)
        except Exception as e:  # noqa: BLE001
            print(f"  警告: {g.locus_tag} をスキップしました ({e})", file=sys.stderr)

    print("[3/5] RNA-seqカバレッジを取得・評価中...", file=sys.stderr)
    coverage_scores_by_gene: Dict[str, Dict[int, float]] = {}
    for gc in all_gc:
        # 候補領域全体をまとめて1回のsamtools depth呼び出しで取得(遺伝子ごとに1回)
        positions = [c.genomic_pos for c in gc.candidates]
        pad = COVERAGE_WINDOW + 5
        region_start = max(1, min(positions) - pad)
        region_end = max(positions) + pad
        try:
            depth = get_pooled_depth(args.bam, gc.gene.seqid, region_start, region_end, args.samtools_bin)
            coverage_scores_by_gene[gc.gene.locus_tag] = score_candidates_by_coverage(gc, depth)
        except subprocess.CalledProcessError as e:
            print(f"  警告: {gc.gene.locus_tag} のカバレッジ取得に失敗 ({e})", file=sys.stderr)

    best_blast_hits: Dict[str, dict] = {}
    if args.orthologs:
        print("[4/5] BLASTPでオルソログと比較中...", file=sys.stderr)
        query_fasta = str(Path(args.workdir) / "extended_proteins.faa")
        blast_tsv = str(Path(args.workdir) / "blast_results.tsv")
        write_extended_fasta(all_gc, query_fasta)
        run_blast(
            query_fasta, args.orthologs, blast_tsv,
            makeblastdb_bin=args.makeblastdb_bin, blastp_bin=args.blastp_bin, threads=args.threads,
        )
        best_blast_hits = parse_best_blast_hits(blast_tsv)
    else:
        print("[4/5] --orthologs 未指定のためBLASTをスキップします", file=sys.stderr)

    print("[5/5] 判定をまとめて出力中...", file=sys.stderr)
    rows = []
    for gc in all_gc:
        cov = coverage_scores_by_gene.get(gc.gene.locus_tag)
        hit = best_blast_hits.get(gc.gene.locus_tag)
        rows.append(make_verdict(gc, cov, hit))

    fieldnames = list(rows[0].keys()) if rows else []
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    n_changed = sum(1 for r in rows if r["verdict"] != "現行アノテーションを支持")
    print(f"完了: {args.out} に {len(rows)} 遺伝子分を出力しました。", file=sys.stderr)
    print(f"  うち現行アノテーションと異なる可能性が示唆された遺伝子: {n_changed}", file=sys.stderr)


if __name__ == "__main__":
    main()

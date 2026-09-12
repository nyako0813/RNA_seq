"""find_true_start_codons.py の中核ロジックに対する簡易テスト(外部ツール非依存)。
samtools / blastp を使わずに済む部分(GFFパース、候補列挙、翻訳、
カバレッジスコアリング、BLAST結果パース、判定ロジック)を検証する。
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from find_true_start_codons import (
    Gene, load_fasta, parse_gff_genes, find_candidates, build_genes_by_seqid_strand,
    translate_cds, get_codon, revcomp, score_candidates_by_coverage,
    parse_best_blast_hits, blast_supported_candidate, make_verdict,
)

def test_revcomp():
    assert revcomp("ATGC") == "GCAT"
    print("OK: revcomp")

def test_translate_cds_alt_start():
    # GTGを開始コドンとして使った場合、最初はMになる
    assert translate_cds("GTGAAATAA") == "MK*"
    assert translate_cds("ATGAAATAA") == "MK*"
    print("OK: translate_cds (alt start -> M)")

def test_find_candidates_plus_strand():
    # 設計: 上流にin-frameのGTGがあり、その手前にin-frameのストップコドンがある
    # 位置(1-based): ...[stop TAA][gap9][GTG candidate][gap9][ATG current]...
    # 読み枠を3の倍数間隔で揃える
    #      1234567890123456789012345678
    genome_seq = "TAAAAACCCGTGAAACCCATGAAATAAGGG"
    #             123456789012345678901234567890
    # 位置10-12 = GTG, 位置19-21 = ATG, 位置1-3 = TAA(上流ストップ)
    assert genome_seq[9:12] == "GTG"
    assert genome_seq[18:21] == "ATG"
    assert genome_seq[0:3] == "TAA"

    genome = {"chr1": genome_seq}
    gene = Gene(locus_tag="geneA", seqid="chr1", start=19, end=27, strand="+")  # ATG..TAA (19-27)
    genes_by_ss = build_genes_by_seqid_strand([gene])

    gc = find_candidates(gene, genome, genes_by_ss, max_upstream=300)
    assert gc.upstream_limit_reason == "stop_codon", gc.upstream_limit_reason
    positions = [c.genomic_pos for c in gc.candidates]
    codons = [c.codon for c in gc.candidates]
    print("candidates:", list(zip(positions, codons)))
    assert 10 in positions, "GTG候補(位置10)が見つかっていない"
    assert 19 in positions, "現行ATG(位置19)が候補に含まれていない"
    # 上流順にソートされているはず(GTGが先、ATGが最後)
    assert positions[0] == 10
    assert positions[-1] == 19
    # 延長タンパク質はGTG(M扱い)から翻訳される
    assert gc.extended_aa_seq.startswith("M")
    print("OK: find_candidates (plus strand)")
    return gc

def test_find_candidates_minus_strand():
    # マイナス鎖の遺伝子: 開始コドン(ATG)はgene.end側(ゲノム座標が大きい方)にあり、
    # 上流(-鎖的に5'方向)はゲノム座標としてはさらに大きい方向に進む。
    # 現行開始コドン位置(pos_1based)=24 のとき、1コドン上流はpos=27、2コドン上流はpos=30。
    # get_codon(strand="-", pos)は seq[pos-3:pos] のrevcompを返すので、
    # 各コドンを置きたい位置の「plus鎖上の実際の塩基」はrevcompして逆算する。
    total_len = 40
    seq_list = ["N"] * total_len  # 1-based用にN埋め、後で0-indexedスライスに使う

    def place_codon(pos_1based: int, codon_on_minus_strand: str):
        plus_strand_triplet = revcomp(codon_on_minus_strand)
        seq_list[pos_1based - 3:pos_1based] = list(plus_strand_triplet)

    place_codon(24, "ATG")  # 現行の開始コドン
    place_codon(27, "GTG")  # 1コドン上流の候補
    place_codon(30, "TAA")  # 2コドン上流のストップ(ここで探索が止まるはず)
    # 残りのN埋め部分をAで埋めて有効な塩基列にしておく(ストップ/候補コドンに偶然一致しないよう配慢)
    genome_seq = "".join(c if c != "N" else "A" for c in seq_list)

    # 念のため意図通りに配置できたか確認
    assert get_codon(genome_seq, genome_seq, 24, "-") == "ATG"
    assert get_codon(genome_seq, genome_seq, 27, "-") == "GTG"
    assert get_codon(genome_seq, genome_seq, 30, "-") == "TAA"

    genome = {"chr1": genome_seq}
    gene = Gene(locus_tag="geneB", seqid="chr1", start=1, end=24, strand="-")
    genes_by_ss = build_genes_by_seqid_strand([gene])
    gc = find_candidates(gene, genome, genes_by_ss, max_upstream=300)
    positions = [c.genomic_pos for c in gc.candidates]
    codons = [c.codon for c in gc.candidates]
    print("minus-strand candidates:", list(zip(positions, codons)))
    assert 24 in positions  # 現行の開始コドン
    assert 27 in positions  # 上流のGTG候補
    assert gc.upstream_limit_reason == "stop_codon"
    assert positions[0] == 27 and positions[-1] == 24  # 上流順ソート確認
    print("OK: find_candidates (minus strand)")

def test_neighbor_gene_boundary():
    # 上流にストップコドンが無くても、隣接遺伝子の終端で探索を打ち切ることを確認
    genome_seq = "ATGAAACCCGTGAAACCCATGAAATAAGGG"
    #             123456789012345678901234567890
    # 位置10-12=GTG, 位置19-21=ATG(現行), 隣接遺伝子が位置1-9に存在すると仮定
    genome = {"chr1": genome_seq}
    gene = Gene(locus_tag="geneA", seqid="chr1", start=19, end=27, strand="+")
    neighbor = Gene(locus_tag="geneUp", seqid="chr1", start=1, end=9, strand="+")
    genes_by_ss = build_genes_by_seqid_strand([gene, neighbor])
    gc = find_candidates(gene, genome, genes_by_ss, max_upstream=300)
    assert gc.upstream_limit_reason == "neighbor_gene", gc.upstream_limit_reason
    positions = [c.genomic_pos for c in gc.candidates]
    assert 10 in positions  # GTGはneighbor(end=9)より後ろなのでまだ候補になる
    print("OK: neighbor gene boundary")

def test_coverage_scoring_and_verdict():
    gc = test_find_candidates_plus_strand()
    # 合成カバレッジ: GTG(位置10)より手前はほぼ0、GTG以降は高い深度
    depth = {}
    for pos in range(1, 40):
        if pos < 10:
            depth[pos] = 1
        else:
            depth[pos] = 50
    scores = score_candidates_by_coverage(gc, depth, window=5)
    print("coverage scores:", scores)
    # GTG(位置10)のスコアが最も高くなるはず
    best_pos = max(scores, key=lambda p: scores[p])
    assert best_pos == 10, best_pos

    # BLASTヒットなしでもRNA-seqだけでverdictが出せることを確認
    verdict = make_verdict(gc, scores, None)
    print("verdict (RNA-seq only):", verdict)
    assert verdict["verdict"].startswith("中確信度") or verdict["verdict"].startswith("高確信度")
    assert verdict["proposed_start_pos"] == 10
    print("OK: coverage scoring + verdict (RNA-seq only)")

    # BLASTヒットもGTGを支持する場合 -> 高確信度で一致するはず
    fake_hit = {"sseqid": "ortholog1", "pident": 95.0, "length": 10, "qstart": 1, "qend": 10, "evalue": 1e-20, "bitscore": 50.0}
    verdict2 = make_verdict(gc, scores, fake_hit)
    print("verdict (RNA-seq + BLAST agree):", verdict2)
    assert verdict2["verdict"] == "高確信度: 上流への延長を提案(RNA-seq・BLAST一致)"
    print("OK: verdict combination (agreement)")

def test_parse_best_blast_hits(tmp_path=None):
    if tmp_path is None:
        tmp_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fake_blast.tsv")
    content = (
        "geneA\tortholog1\t90.0\t50\t1\t50\t1\t50\t1e-30\t100.0\n"
        "geneA\tortholog2\t80.0\t50\t3\t50\t1\t48\t1e-10\t60.0\n"
        "geneC\tortholog3\t99.0\t20\t1\t20\t1\t20\t1e-15\t80.0\n"
    )
    with open(tmp_path, "w") as fh:
        fh.write(content)
    best = parse_best_blast_hits(tmp_path)
    assert best["geneA"]["sseqid"] == "ortholog1"  # bitscore高い方
    assert best["geneA"]["qstart"] == 1
    assert best["geneC"]["sseqid"] == "ortholog3"
    print("OK: parse_best_blast_hits")


if __name__ == "__main__":
    test_revcomp()
    test_translate_cds_alt_start()
    test_find_candidates_plus_strand()
    test_find_candidates_minus_strand()
    test_neighbor_gene_boundary()
    test_coverage_scoring_and_verdict()
    test_parse_best_blast_hits()
    print("\nすべてのテストに合格しました。")

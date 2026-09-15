"""resume_gca_gcf_step9e.py の中核ロジック(BLAST/samtools呼び出しを伴わない部分)に対する簡易テスト。"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gca_gcf_comparison import GeneRecord
from resume_gca_gcf_step9e import (
    recommend_frame,
    orf_span_from_gca_start,
    check_gca_start_validity,
    step_ratio_at,
)


def hit(pident, qstart):
    return {"pident": pident, "qstart": qstart}


def test_recommend_frame_only_gcf_good():
    assert recommend_frame(hit(90.0, 1), None) == "GCF読み枠が妥当"
    print("OK: recommend_frame (GCFのみヒットあり)")


def test_recommend_frame_only_gca_good():
    assert recommend_frame(None, hit(90.0, 1)) == "GCA読み枠が妥当"
    print("OK: recommend_frame (GCAのみヒットあり)")


def test_recommend_frame_both_weak_or_none():
    assert recommend_frame(None, None) == "判断つかない(両方ヒットが弱い、または拮抗)"
    assert recommend_frame(hit(10.0, 50), hit(12.0, 60)) == "判断つかない(両方ヒットが弱い、または拮抗)"
    print("OK: recommend_frame (ヒットなし/両方弱い)")


def test_recommend_frame_both_good_prefers_higher_pident():
    assert recommend_frame(hit(95.0, 1), hit(60.0, 1)) == "GCF読み枠が妥当"
    assert recommend_frame(hit(60.0, 1), hit(95.0, 1)) == "GCA読み枠が妥当"
    print("OK: recommend_frame (両方良好、pident差で決着)")


def test_recommend_frame_both_good_tie():
    assert recommend_frame(hit(80.0, 1), hit(82.0, 1)) == "判断つかない(両方ヒットが弱い、または拮抗)"
    print("OK: recommend_frame (両方良好、拮抗)")


def test_orf_span_from_gca_start_plus_strand():
    gene = GeneRecord(locus_tag="g1", start=100, end=300, strand="+")
    assert orf_span_from_gca_start(gene, 91) == (91, 300)
    print("OK: orf_span_from_gca_start (+鎖)")


def test_orf_span_from_gca_start_minus_strand():
    gene = GeneRecord(locus_tag="g1", start=100, end=300, strand="-")
    assert orf_span_from_gca_start(gene, 309) == (100, 309)
    print("OK: orf_span_from_gca_start (-鎖)")


def test_check_gca_start_validity_valid_plus_strand():
    # 9bp上流にATG、間にストップコドンなし、現行のストップ(TAA)まで妥当なORF
    seq = "N" * 9 + "ATG" + "AAA" + "TAA" + "N" * 20
    # 候補開始(lo)=10 (1-based, "ATG"の"A"), 現行ストップまで(hi)=18 ("TAA"の"A")
    valid, reason = check_gca_start_validity(seq, lo=10, hi=18, strand="+")
    assert valid, reason
    print("OK: check_gca_start_validity (+鎖, 妥当)")


def test_check_gca_start_validity_internal_stop():
    seq = "N" * 9 + "ATG" + "TAA" + "ATG" + "TAA" + "N" * 10
    valid, reason = check_gca_start_validity(seq, lo=10, hi=21, strand="+")
    assert not valid
    assert "ストップコドン" in reason
    print("OK: check_gca_start_validity (+鎖, 読み枠内ストップコドンで無効)")


def test_check_gca_start_validity_noncanonical_codon():
    seq = "N" * 9 + "CCC" + "AAA" + "TAA" + "N" * 10
    valid, reason = check_gca_start_validity(seq, lo=10, hi=18, strand="+")
    assert not valid
    assert "非正準" in reason
    print("OK: check_gca_start_validity (+鎖, 非正準コドンで無効)")


def test_step_ratio_at_plus_strand_detects_rise():
    # posの直前30bpはdepth 0、直後30bpはdepth 100 -> ステップ比が大きいはず
    depth = {}
    pos = 1000
    for p in range(pos - 30, pos):
        depth[p] = 0
    for p in range(pos, pos + 30):
        depth[p] = 100
    ratio = step_ratio_at(depth, pos, "+")
    assert ratio > 50
    print("OK: step_ratio_at (+鎖, 立ち上がり検出)")


def test_step_ratio_at_minus_strand_detects_rise():
    # -鎖では「直前(下流側, 座標が大きい方)」に発現がないことがカバレッジの立ち上がりに相当
    depth = {}
    pos = 1000
    for p in range(pos + 1, pos + 31):
        depth[p] = 0
    for p in range(pos - 29, pos + 1):
        depth[p] = 100
    ratio = step_ratio_at(depth, pos, "-")
    assert ratio > 50
    print("OK: step_ratio_at (-鎖, 立ち上がり検出)")


if __name__ == "__main__":
    test_recommend_frame_only_gcf_good()
    test_recommend_frame_only_gca_good()
    test_recommend_frame_both_weak_or_none()
    test_recommend_frame_both_good_prefers_higher_pident()
    test_recommend_frame_both_good_tie()
    test_orf_span_from_gca_start_plus_strand()
    test_orf_span_from_gca_start_minus_strand()
    test_check_gca_start_validity_valid_plus_strand()
    test_check_gca_start_validity_internal_stop()
    test_check_gca_start_validity_noncanonical_codon()
    test_step_ratio_at_plus_strand_detects_rise()
    test_step_ratio_at_minus_strand_detects_rise()
    print("\nすべてのテストに合格しました。")

"""gca_gcf_comparison.py / integrate_gca_evidence.py の中核ロジックに対する簡易テスト。"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gca_gcf_comparison import GeneRecord, translation_start, build_comparison, parse_old_locus_tag
from integrate_gca_evidence import classify_gca_agreement, check_upstream_extension_valid, translate


def test_translation_start_plus_minus():
    plus = GeneRecord(locus_tag="a", start=100, end=200, strand="+")
    minus = GeneRecord(locus_tag="b", start=100, end=200, strand="-")
    assert translation_start(plus) == 100
    assert translation_start(minus) == 200
    print("OK: translation_start")


def test_parse_old_locus_tag_picks_underscore_form():
    assert parse_old_locus_tag("MA0001%2CMA_0001") == "MA_0001"
    assert parse_old_locus_tag(None) is None
    print("OK: parse_old_locus_tag")


def test_build_comparison_diff_direction_plus_strand():
    # +鎖: GCAがGCFより上流(小さい座標) -> diff_bpは正(延長方向)
    gcf = {"g1": GeneRecord("g1", start=100, end=300, strand="+", old_locus_tag="old1")}
    gca = {"old1": GeneRecord("old1", start=91, end=300, strand="+")}
    rows = build_comparison(gcf, gca)
    assert rows[0]["diff_bp"] == 9
    assert rows[0]["diff_codons"] == 3.0
    assert rows[0]["frame_shift_flag"] is False
    print("OK: build_comparison (+鎖, 延長方向)")


def test_build_comparison_diff_direction_minus_strand():
    # -鎖: GCAがGCFより上流(大きい座標) -> diff_bpは正(延長方向)
    gcf = {"g1": GeneRecord("g1", start=100, end=300, strand="-", old_locus_tag="old1")}
    gca = {"old1": GeneRecord("old1", start=100, end=309, strand="-")}
    rows = build_comparison(gcf, gca)
    assert rows[0]["diff_bp"] == 9
    print("OK: build_comparison (-鎖, 延長方向)")


def test_build_comparison_frame_shift_flag():
    gcf = {"g1": GeneRecord("g1", start=100, end=300, strand="+", old_locus_tag="old1")}
    gca = {"old1": GeneRecord("old1", start=95, end=300, strand="+")}
    rows = build_comparison(gcf, gca)
    assert rows[0]["frame_shift_flag"] is True
    print("OK: build_comparison (frame_shift_flag)")


def test_classify_gca_agreement():
    assert classify_gca_agreement(100, 90, 90, None) == "GCA比較不可"
    assert classify_gca_agreement(100, 90, 90, 100) == "現行(GCF)と一致"
    assert classify_gca_agreement(100, 90, 80, 90) == "RNA-seqと一致"
    assert classify_gca_agreement(100, 90, 80, 80) == "BLASTと一致"
    assert classify_gca_agreement(100, 90, 90, 90) == "BLAST・RNA-seq両方と一致"
    assert classify_gca_agreement(100, 90, 80, 70) == "いずれとも異なる(第3の候補)"
    print("OK: classify_gca_agreement")


def test_check_upstream_extension_valid():
    # +鎖: ATG(候補)...ATG(現行)で間にストップコドンなし
    seq = "N" * 9 + "ATG" + "AAA" + "ATG" + "N" * 20
    # 候補位置=10(1-based, ATG開始), 現行位置=16(ATG)
    valid, reason = check_upstream_extension_valid(seq, current_start=16, candidate_start=10, strand="+")
    assert valid, reason

    # 間にストップコドン(TAA)がある場合は無効
    seq2 = "N" * 9 + "ATG" + "TAA" + "ATG" + "N" * 20
    valid2, reason2 = check_upstream_extension_valid(seq2, current_start=16, candidate_start=10, strand="+")
    assert not valid2

    # 候補位置自体が非正準コドンなら無効
    seq3 = "N" * 9 + "CCC" + "AAA" + "ATG" + "N" * 20
    valid3, reason3 = check_upstream_extension_valid(seq3, current_start=16, candidate_start=10, strand="+")
    assert not valid3
    print("OK: check_upstream_extension_valid")


def test_translate_stop_codon():
    assert translate("ATGAAATAA") == "MK*"
    print("OK: translate")


if __name__ == "__main__":
    test_translation_start_plus_minus()
    test_parse_old_locus_tag_picks_underscore_form()
    test_build_comparison_diff_direction_plus_strand()
    test_build_comparison_diff_direction_minus_strand()
    test_build_comparison_frame_shift_flag()
    test_classify_gca_agreement()
    test_check_upstream_extension_valid()
    test_translate_stop_codon()
    print("\nすべてのテストに合格しました。")

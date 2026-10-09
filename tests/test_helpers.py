import pandas as pd

from tast.filters import choose_strand, merge_intervals, transcript_stats
from tast.layout import colinear_runs, pack_rows, region_labels, region_order, spread_positions
from tast.orfs import find_orfs


def test_find_orfs_basic():
    #ATG AAA TAG = 9 nt ORF at 1-based 3..11
    seq = "CC" + "ATGAAATAG" + "GG"
    assert find_orfs(seq, min_len=0) == ([[3, 11]], [3, 11])
    assert find_orfs(seq, min_len=9) == ([], None)


def test_find_orfs_internal_atg_and_no_stop():
    seq = "ATGATGTTTTAA" + "CATGCCC"  #second frame-0 ATG is internal; last ATG has no stop
    all_orfs, longest = find_orfs(seq, min_len=0)
    assert all_orfs == [[1, 12], [4, 12]]
    assert longest == [1, 12]


def test_find_orfs_lowercase():
    assert find_orfs("ccatgaaatagg", 0) == find_orfs("CCATGAAATAGG", 0)


def test_region_labels_transcript_position():
    assert region_labels([10, 300, 500]) == [0, 1, 2]          #forward strand gene
    assert region_labels([500, 300, 10]) == [2, 1, 0]          #reverse strand gene
    assert region_labels([890, 170]) == [1, 0]
    assert region_labels([500, 10, 500, 20]) == [2, 0, 2, 1]   #repeated part shares a number


def test_colinear_runs_and_order():
    assert colinear_runs([0, 1, 2, 3]) == [(0, 3, 1)]
    assert colinear_runs([3, 2, 1, 0]) == [(0, 3, -1)]
    assert colinear_runs([0, 1, 4, 5, 2, 3]) == [(0, 3, 1), (4, 5, 1)]   #block moved
    assert colinear_runs([0, 1, 2, 2, 3]) == [(0, 2, 1), (3, 4, 1)]      #repeat breaks the run
    assert colinear_runs([0]) == [(0, 0, 0)]
    assert region_order([0, 1, 2]) == "forward"
    assert region_order([2, 1, 0]) == "reverse"
    assert region_order([0, 1, 4, 5, 2, 3]) == "mixed"
    assert region_order([0]) == "single"


def test_spread_positions_does_not_modify_input():
    pos = [0, 100, 1000]
    assert spread_positions(pos) == [0, 250, 1000]
    assert pos == [0, 100, 1000]


def test_pack_rows_uses_minimum_rows():
    starts = [0, 5, 20, 25, 40]
    ends = [10, 30, 35, 50, 45]
    rows = pack_rows(starts, ends)
    #the minimum possible number of rows is the largest number of transcripts overlapping one position
    depth = max(sum(s <= x <= e for s, e in zip(starts, ends)) for x in range(0, 51))
    assert depth == 3
    assert len(rows) == depth
    for row in rows:
        for a, b in zip(row, row[1:]):
            assert ends[a] < starts[b]


def test_pack_rows_min_gap():
    assert len(pack_rows([0, 15], [10, 20], min_gap=0)) == 1
    assert len(pack_rows([0, 15], [10, 20], min_gap=10)) == 2


def _blast(rows):
    cols = ["qseqid", "qlen", "qstart", "qend", "sstart", "send", "strand", "bitscore"]
    return pd.DataFrame(rows, columns=cols)


def test_choose_strand_dominant():
    b = _blast([["t1", 100, 1, 50, 10, 60, "+", 80.0],
                ["t1", 100, 51, 100, 200, 250, "-", 20.0],
                ["t2", 100, 1, 50, 10, 60, "-", 50.0]])
    kept = choose_strand(b, "dominant")
    assert kept[["qseqid", "strand"]].values.tolist() == [["t1", "+"], ["t2", "-"]]
    assert len(choose_strand(b, "both")) == 3


def test_transcript_stats_coverage_and_intron():
    b = _blast([["t1", 100, 1, 40, 100, 140, "+", 1.0],
                ["t1", 100, 31, 60, 1000, 1030, "+", 1.0],   #overlaps on transcript, 860 bp gap on genome
                ["t2", 200, 1, 50, 5, 55, "+", 1.0]])
    s = transcript_stats(b)
    assert s.loc["t1", "coverage"] == 59.0       #(40-1) + (60-40) merged = 59 of 100
    assert s.loc["t1", "max_intron"] == 860
    assert s.loc["t2", "max_intron"] == 0
    assert s.loc["t2", "coverage"] == 24.5


def test_merge_intervals_touching_merge():
    df = pd.DataFrame({"qseqid": ["a"] * 3, "s": [1, 10, 30], "e": [10, 20, 40]})
    m = merge_intervals(df, "s", "e")
    assert m[["Left", "Right"]].values.tolist() == [[1, 20], [30, 40]]

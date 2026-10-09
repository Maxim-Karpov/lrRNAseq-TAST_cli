"""Strand selection, transcript coverage / intron statistics and filtering."""
from __future__ import annotations

import pandas as pd

STRAND_MODES = ("dominant", "both", "plus", "minus")


def choose_strand(blast: pd.DataFrame, mode: str) -> pd.DataFrame:
    """Keep alignments on the requested strand(s).

    dominant - for each transcript, the strand with the highest total bitscore (ties go to "+")
    both     - all alignments (behaviour of the original notebook, version 1.0)
    plus     - only alignments hitting the genome in forward orientation
    minus    - only alignments hitting the genome in reverse orientation
    """
    if mode == "both":
        kept = blast
    elif mode == "plus":
        kept = blast[blast["strand"] == "+"]
    elif mode == "minus":
        kept = blast[blast["strand"] == "-"]
    elif mode == "dominant":
        score = blast.groupby(["qseqid", "strand"])["bitscore"].sum().unstack(fill_value=0)
        for s in ("+", "-"):
            if s not in score:
                score[s] = 0
        dominant = (score["-"] > score["+"]).map({True: "-", False: "+"})
        kept = blast[blast["strand"] == blast["qseqid"].map(dominant)]
    else:
        raise ValueError(f"strand mode must be one of {STRAND_MODES}")
    return kept.reset_index(drop=True)


def merge_intervals(df: pd.DataFrame, start_col: str, end_col: str,
                    group_col: str = "qseqid") -> pd.DataFrame:
    """For each group, merge overlapping or touching [start, end] intervals.

    Returns one row per merged block with columns group_col, "block", "Left", "Right".
    """
    df = df.sort_values([group_col, start_col])
    running_end = df.groupby(group_col)[end_col].cummax()
    previous_end = running_end.groupby(df[group_col]).shift()
    new_block = (df[start_col] > previous_end) | previous_end.isna()
    block = new_block.astype(int).groupby(df[group_col]).cumsum()
    merged = (df.assign(block=block)
                .groupby([group_col, "block"])
                .agg(Left=(start_col, "min"), Right=(end_col, "max")))
    return merged.reset_index()


def transcript_stats(blast: pd.DataFrame) -> pd.DataFrame:
    """Per-transcript statistics, indexed by qseqid.

    coverage    - % of the transcript length covered by the union of its alignments
    max_intron  - largest gap on the genome between merged alignment blocks (0 if one block)
    """
    merged_q = merge_intervals(blast, "qstart", "qend")
    covered = (merged_q["Right"] - merged_q["Left"]).groupby(merged_q["qseqid"]).sum()
    grouped = blast.groupby("qseqid")
    qlen = grouped["qlen"].max()

    merged_s = merge_intervals(blast, "sstart", "send")
    merged_s["gap"] = merged_s["Left"] - merged_s.groupby("qseqid")["Right"].shift()
    max_gap = merged_s.groupby("qseqid")["gap"].max().fillna(0).astype(int)

    stats = pd.DataFrame({
        "qlen": qlen,
        "n_alignments": grouped.size(),
        "strand": grouped["strand"].agg(lambda s: "".join(sorted(set(s), reverse=True))),
        "coverage": 100 / qlen * covered,
        "max_intron": max_gap,
    })
    stats.index.name = "qseqid"
    return stats


def apply_filters(blast: pd.DataFrame, stats: pd.DataFrame, max_intron: int,
                  min_coverage: float) -> pd.DataFrame:
    """Return the alignments of transcripts with max_intron < max_intron and
    coverage > min_coverage, sorted along the genome."""
    passed = stats.index[(stats["max_intron"] < max_intron) & (stats["coverage"] > min_coverage)]
    kept = blast[blast["qseqid"].isin(set(passed))]
    return kept.sort_values(["sstart", "qstart"], kind="stable").reset_index(drop=True)

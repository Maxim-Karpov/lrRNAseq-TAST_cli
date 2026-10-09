"""The TAST analysis: from BLAST alignments to plot rows."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional

import pandas as pd

from .filters import apply_filters, choose_strand, merge_intervals, transcript_stats
from .files import InputError, subset_fasta
from .layout import colinear_runs, pack_rows, region_labels, region_order, spread_positions
from .orfs import find_orfs


@dataclass
class Settings:
    max_intron: int = 20000
    min_coverage: float = 25.0
    strand: str = "dominant"
    barcode_len: int = 0
    min_orf_len: int = 26
    min_gap: int = 0


@dataclass
class Row:
    """One row of the plot."""
    transcripts: List[str]          #transcript names, left to right
    ranges: List[List[int]]         #genomic [start, end] of each transcript
    regions: List[List[int]]        #genomic [start, end] of every alignment region in the row
    in_orf: List[bool]              #region overlaps any ORF of its transcript
    in_largest_orf: List[bool]      #region overlaps the largest ORF of its transcript
    labels: List[int]               #region numbers drawn on the plot
    label_x: List[int]              #x positions of the numbers
    frames: List[List[int]]         #genomic [start, end] of colinear blocks of regions
    frame_directions: List[int]     #1 forward, -1 reverse, 0 single region


@dataclass
class Result:
    rows: List[Row]
    summary: pd.DataFrame           #every transcript in the BLAST file with its statistics and status
    filtered_ids: List[str]
    warnings: List[str] = field(default_factory=list)


def _overlaps(a: List[int], b: List[int]) -> bool:
    """True if the half-open ranges [a0, a1) and [b0, b1) share a position."""
    return max(a[0], b[0]) < min(a[1], b[1])


def run(blast: pd.DataFrame, transcripts_fasta: Path, filtered_fasta_out: Path,
        settings: Settings, log: Callable[[str], None] = print) -> Result:
    warnings: List[str] = []

    def warn(msg: str) -> None:
        warnings.append(msg)
        log("WARNING: " + msg)

    #---- strand choice, coverage and intron filters
    summary = pd.DataFrame(index=pd.Index(sorted(blast["qseqid"].unique()), name="qseqid"))
    stranded = choose_strand(blast, settings.strand)
    stats = transcript_stats(stranded)
    summary = summary.join(stats)
    summary["status"] = "removed: no alignments on chosen strand"
    has_stats = summary["coverage"].notna()
    summary.loc[has_stats & (summary["max_intron"] >= settings.max_intron), "status"] = \
        f"removed: max intron >= {settings.max_intron}"
    summary.loc[has_stats & (summary["coverage"] <= settings.min_coverage), "status"] = \
        f"removed: coverage <= {settings.min_coverage:g}%"

    kept = apply_filters(stranded, stats, settings.max_intron, settings.min_coverage)
    filtered_ids = kept.drop_duplicates("qseqid")["qseqid"].tolist()
    log(f"{len(filtered_ids)} of {summary.shape[0]} transcripts pass the strand, intron and coverage filters")
    if not filtered_ids:
        raise InputError("No transcripts passed the filters. Try a larger --max-intron or a smaller --min-coverage.")

    #---- sequences of the filtered transcripts
    seqs = subset_fasta(transcripts_fasta, set(filtered_ids), filtered_fasta_out)
    missing = [i for i in filtered_ids if i not in seqs]
    if missing:
        warn(f"{len(missing)} filtered transcripts were not found in {transcripts_fasta} "
             f"(e.g. {missing[0]}); they are not plotted")
        summary.loc[missing, "status"] = "removed: not in transcripts FASTA"

    #---- ORFs (found on barcode-trimmed sequences, reported in BLAST coordinates)
    b = settings.barcode_len
    qlen = kept.groupby("qseqid")["qlen"].max()
    orfs: Dict[str, List[List[int]]] = {}
    largest: Dict[str, List[int]] = {}
    odd_offsets = []
    for tid, seq in seqs.items():
        trimmed = seq[b:-b] if b else seq
        all_orfs, longest = find_orfs(trimmed, settings.min_orf_len)
        if longest is None:
            summary.loc[tid, "status"] = "removed: no ORF"
            continue
        offset = (int(qlen[tid]) - len(trimmed)) // 2
        if offset not in (0, b):
            odd_offsets.append(tid)
        orfs[tid] = [[c + offset for c in orf] for orf in all_orfs]
        largest[tid] = [c + offset for c in longest]
    no_orf = [t for t in seqs if t not in orfs]
    if no_orf:
        log(f"{len(no_orf)} transcripts have no ORF longer than {settings.min_orf_len} nt and are not plotted")
    if odd_offsets:
        warn(f"for {len(odd_offsets)} transcripts the length in the BLAST file differs from the FASTA length "
             f"by something other than twice the barcode length (e.g. {odd_offsets[0]}); "
             "check --barcode-len and that BLAST was run on these sequences")
    if not orfs:
        raise InputError("No ORFs were found in any filtered transcript.")

    #---- one model per transcript: alignment regions in genomic order
    models = []
    for tid, group in kept.groupby("qseqid", sort=True):
        if tid not in orfs:
            continue
        regions = [[s, e] for s, e in zip(group["sstart"].tolist(), group["send"].tolist())]
        qranges = [[s, e] for s, e in zip(group["qstart"].tolist(), group["qend"].tolist())]
        orf_blocks = merge_intervals(pd.DataFrame({"qseqid": tid, "Left": [o[0] for o in orfs[tid]],
                                                   "Right": [o[1] for o in orfs[tid]]}), "Left", "Right")
        orf_blocks = orf_blocks[["Left", "Right"]].values.tolist()
        start = min(c for r in regions for c in r)
        end = max(c for r in regions for c in r)
        models.append({
            "qseqid": tid, "start": start, "end": end, "TrLen": end - start, "regions": regions,
            "in_orf": [any(_overlaps(q, o) for o in orf_blocks) for q in qranges],
            "in_largest_orf": [_overlaps(q, largest[tid]) for q in qranges],
            "labels": region_labels([q[0] for q in qranges]),
        })
    models = pd.DataFrame(models).sort_values(by=["start", "end"]).reset_index(drop=True)

    #---- rows
    members = pack_rows(models["start"].tolist(), models["end"].tolist(), settings.min_gap)
    rows = pd.DataFrame({
        "members": members,
        "Total_aln_len": [int(models["TrLen"][m].sum()) for m in members],
        "span": [int(models["end"][m].max() - models["start"][m].min()) for m in members],
        "first_start": [int(models["start"][m].min()) for m in members],
    })
    #widest rows first (as in the original notebook); ties broken explicitly so the order is
    #the same on every computer and library version
    rows = rows.sort_values(["span", "Total_aln_len", "first_start"], ascending=[False, False, True],
                            kind="stable").reset_index(drop=True)

    plot_rows: List[Row] = []
    summary["plot_row"] = pd.NA
    for row_number, m in enumerate(rows["members"]):
        sub = models.iloc[m]
        regions = [r for rs in sub["regions"] for r in rs]
        frames, directions, offset = [], [], 0
        for rs, labels in zip(sub["regions"], sub["labels"]):
            for first, last, direction in colinear_runs(labels):
                frames.append([regions[first + offset][0], regions[last + offset][1]])
                directions.append(direction)
            offset += len(rs)
        plot_rows.append(Row(
            transcripts=sub["qseqid"].tolist(),
            ranges=[[s, e] for s, e in zip(sub["start"], sub["end"])],
            regions=regions,
            in_orf=[f for fs in sub["in_orf"] for f in fs],
            in_largest_orf=[f for fs in sub["in_largest_orf"] for f in fs],
            labels=[x for ls in sub["labels"] for x in ls],
            label_x=spread_positions([int(r[0] + ((r[1] - r[0]) / 2)) for r in regions]),
            frames=frames,
            frame_directions=directions,
        ))
        summary.loc[sub["qseqid"].tolist(), ["status", "plot_row"]] = ["plotted", row_number]

    #largest ORF in BLAST transcript coordinates, for the output table
    summary["largest_orf"] = pd.Series({t: f"{c[0]}-{c[1]}" for t, c in largest.items()}, dtype=object)
    summary["n_orfs"] = pd.Series({t: len(o) for t, o in orfs.items()}, dtype="Int64")
    summary["region_numbers"] = pd.Series(
        {t: ",".join(map(str, ls)) for t, ls in zip(models["qseqid"], models["labels"])}, dtype=object)
    summary["region_order"] = pd.Series(
        {t: region_order(ls) for t, ls in zip(models["qseqid"], models["labels"])}, dtype=object)
    return Result(rows=plot_rows, summary=summary, filtered_ids=filtered_ids, warnings=warnings)

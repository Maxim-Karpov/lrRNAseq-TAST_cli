"""Arranging transcripts into plot rows and numbering their alignment regions."""
from __future__ import annotations

import heapq
from typing import List, Sequence, Tuple


def pack_rows(starts: Sequence[int], ends: Sequence[int], min_gap: int = 0) -> List[List[int]]:
    """Greedy interval partitioning: put transcripts that do not overlap on the genome in
    the same row, using the smallest possible number of rows.

    Transcripts must be given sorted by start. Each is placed in the row whose last
    transcript ended earliest, if that is more than `min_gap` bp before its start.
    Returns a list of rows, each a list of transcript indices.
    """
    rows: List[List[int]] = []
    free = []  #heap of (end of last transcript in row, row number)
    for i, (start, end) in enumerate(zip(starts, ends)):
        if free and free[0][0] + min_gap < start:
            _, row = heapq.heappop(free)
            rows[row].append(i)
        else:
            row = len(rows)
            rows.append([i])
        heapq.heappush(free, (end, row))
    return rows


def region_labels(qstarts: Sequence[int]) -> List[int]:
    """Number the alignment regions of one transcript by their position in the transcript.
    `qstarts` are the transcript start positions of the regions, given in genomic order.

    0 is the region nearest the transcript's 5' end, so along the genome the numbers read
    0, 1, 2 ... for a gene on the forward strand, ... 2, 1, 0 for a gene on the reverse strand,
    and out of order where parts of the transcript are rearranged. Regions with the same qstart
    (the same part of the transcript aligned in several places) share a number.
    """
    rank = {q: i for i, q in enumerate(sorted(set(qstarts)))}
    return [rank[q] for q in qstarts]


def colinear_runs(labels: Sequence[int]) -> List[Tuple[int, int, int]]:
    """Split region numbers (in genomic order) into colinear blocks: runs where the numbers
    keep increasing (direction 1, transcript in the same orientation as the genome) or keep
    decreasing (direction -1, reverse orientation). A repeated number ends a run. Runs are
    found greedily from the left. Returns (first index, last index, direction) for each run;
    a run of one region has direction 0.
    """
    runs = []
    i, n = 0, len(labels)
    while i < n:
        j, direction = i, 0
        if i + 1 < n and labels[i + 1] != labels[i]:
            direction = 1 if labels[i + 1] > labels[i] else -1
            j = i + 1
            while j + 1 < n and (labels[j + 1] - labels[j]) * direction > 0:
                j += 1
        runs.append((i, j, direction))
        i = j + 1
    return runs


def region_order(labels: Sequence[int]) -> str:
    """Summarise the order of a transcript's regions along the genome:
    forward / reverse (one colinear block), mixed (several blocks) or single (one region)."""
    if len(labels) == 1:
        return "single"
    runs = colinear_runs(labels)
    if len(runs) == 1:
        return "forward" if runs[0][2] == 1 else "reverse"
    return "mixed"


def spread_positions(positions: Sequence[int], spacing: int = 250) -> List[int]:
    """Push label positions apart so neighbours are at least `spacing` bp apart."""
    out = list(positions)
    for i in range(len(out) - 1):
        if out[i + 1] - out[i] < spacing:
            out[i + 1] = out[i] + spacing
    return out

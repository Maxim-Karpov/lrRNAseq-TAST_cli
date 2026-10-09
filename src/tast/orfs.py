"""Open reading frame finding on the forward strand of transcripts."""
from __future__ import annotations

from typing import List, Optional, Tuple

STOPS = ("TAA", "TGA", "TAG")


def find_orfs(seq: str, min_len: int = 26) -> Tuple[List[List[int]], Optional[List[int]]]:
    """Find ORFs (ATG ... first in-frame stop, stop included) in the three forward frames.

    Coordinates are 1-based and inclusive: [start, stop]. Only ORFs longer than
    `min_len` nucleotides are returned. Every ATG gives an ORF, so ORFs starting at
    internal in-frame ATGs are included. ORFs without a stop codon are ignored.

    Returns (all_orfs, longest_orf); longest_orf is None if no ORF passes min_len.
    """
    seq = seq.upper()
    n = len(seq)
    orfs: List[List[int]] = []
    for frame in range(3):
        last = n - (n - frame) % 3  #end of the last complete codon in this frame
        codon_starts = range(frame, last, 3)
        #scan backwards to know, for every codon, where the next in-frame stop ends
        next_stop_end = {}
        stop_end = None
        for i in reversed(codon_starts):
            if seq[i:i + 3] in STOPS:
                stop_end = i + 3
            next_stop_end[i] = stop_end
        for i in codon_starts:
            if seq[i:i + 3] == "ATG" and next_stop_end[i] is not None:
                orfs.append([i + 1, next_stop_end[i]])
    lengths = [stop - start + 1 for start, stop in orfs]
    keep = [k for k, length in enumerate(lengths) if length > min_len]
    if not keep:
        return [], None
    longest_len = max(lengths)
    longest = next(orfs[k] for k in keep if lengths[k] == longest_len)
    return [orfs[k] for k in keep], longest

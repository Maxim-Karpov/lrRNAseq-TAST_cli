"""Reading and writing BLAST tables and FASTA files."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Set, Tuple

import pandas as pd

BLAST_COLUMNS = [
    "qseqid", "sseqid", "slen", "qlen", "qstart", "qend", "sstart", "send",
    "length", "mismatch", "gapopen", "pident", "evalue", "bitscore",
]
#the -outfmt string BLAST must be run with
BLAST_OUTFMT = "6 " + " ".join(BLAST_COLUMNS)


class InputError(Exception):
    """A problem with the user's input files, reported without a traceback."""


def read_blast(path: Path) -> pd.DataFrame:
    """Read a BLAST -outfmt '6 qseqid sseqid slen qlen qstart qend sstart send length
    mismatch gapopen pident evalue bitscore' table.

    Alignments where the genome is hit in reverse orientation (sstart > send) get
    strand "-" and their genomic coordinates are swapped so that sstart < send.
    """
    path = Path(path)
    if not path.is_file():
        raise InputError(f"BLAST file not found: {path}")
    with open(path) as f:
        first = f.readline()
    if not first.strip():
        raise InputError(f"BLAST file is empty: {path}")
    fields = first.rstrip("\r\n").split("\t")
    if len(fields) != len(BLAST_COLUMNS):
        raise InputError(
            f"{path} has {len(fields)} tab-separated columns but {len(BLAST_COLUMNS)} are needed.\n"
            f"Run BLAST with: -outfmt '{BLAST_OUTFMT}'")
    has_header = fields[0] == "qseqid"
    blast = pd.read_csv(path, sep="\t", header=None, names=BLAST_COLUMNS,
                        skiprows=1 if has_header else 0,
                        dtype={"qseqid": str, "sseqid": str})
    for col in ["slen", "qlen", "qstart", "qend", "sstart", "send"]:
        if not pd.api.types.is_integer_dtype(blast[col]):
            raise InputError(f"Column '{col}' of {path} contains non-integer values. "
                             f"Run BLAST with: -outfmt '{BLAST_OUTFMT}'")
    blast["strand"] = "+"
    minus = blast["sstart"] > blast["send"]
    blast.loc[minus, "strand"] = "-"
    blast.loc[minus, ["sstart", "send"]] = blast.loc[minus, ["send", "sstart"]].to_numpy()
    return blast


def iter_fasta(path: Path) -> Iterator[Tuple[str, str, str]]:
    """Yield (id, full header line without '>', sequence) for each FASTA record."""
    path = Path(path)
    if not path.is_file():
        raise InputError(f"FASTA file not found: {path}")
    header = None
    chunks: List[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header.split()[0] if header.split() else "", header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            elif line.strip():
                if header is None:
                    raise InputError(f"{path} does not look like a FASTA file (no '>' header before sequence)")
                chunks.append(line.strip())
    if header is not None:
        yield header.split()[0] if header.split() else "", header, "".join(chunks)


def read_genome(path: Path) -> Tuple[str, str]:
    """Read a FASTA file that must contain exactly one sequence. Returns (id, sequence)."""
    records = list(iter_fasta(path))
    if len(records) == 0:
        raise InputError(f"No sequence found in genome FASTA {path}")
    if len(records) > 1:
        raise InputError(f"Genome FASTA {path} contains {len(records)} sequences; "
                         "TAST works on one genomic region at a time, so give a FASTA with a single sequence")
    return records[0][0], records[0][2]


def subset_fasta(path: Path, keep: Set[str], out_path: Path) -> Dict[str, str]:
    """Write the records of `path` whose id is in `keep` to `out_path` (in input order,
    60 bases per line) and return them as {id: sequence}."""
    found: Dict[str, str] = {}
    with open(out_path, "w", newline="\n") as out:
        for rid, header, seq in iter_fasta(path):
            if rid in keep and rid not in found:
                found[rid] = seq
                out.write(f">{header}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i + 60] + "\n")
    return found


def write_lines(lines: Iterable[str], path: Path) -> None:
    with open(path, "w", newline="\n") as f:
        for line in lines:
            f.write(f"{line}\n")

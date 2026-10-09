"""Command line interface: tast -t transcripts.fasta -g genome.fasta -o results"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import List, Optional

from . import __version__
from .files import BLAST_OUTFMT, InputError, read_blast, read_genome, write_lines
from .filters import STRAND_MODES
from .pipeline import Settings, run

EXAMPLES = """examples:
  run BLAST and plot (needs BLAST+ installed):
    tast -t transcripts.fasta -g region.fasta -o results

  re-plot from an existing BLAST result with stricter filters (fast, no BLAST):
    tast -t transcripts.fasta -b results/blast_hits.tsv -o results_cov50 --min-coverage 50

  try it on the example data from the repository:
    tast -t example/all_transcripts.fasta -g example/Ef318A10.fst -b example/lrRNAseq_to_Ef318A10.out -o example_results
"""


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="tast",
        description="lrRNAseq TAST: align full-length long-read RNAseq transcripts to a genomic region "
                    "(e.g. a ~100 kb BAC) and plot the total alignment to reveal genes and their splicing.",
        epilog=EXAMPLES, formatter_class=argparse.RawDescriptionHelpFormatter)

    inp = p.add_argument_group("input")
    inp.add_argument("-t", "--transcripts", required=True, type=Path, metavar="FASTA",
                     help="long-read transcripts (FASTA)")
    inp.add_argument("-g", "--genome", type=Path, metavar="FASTA",
                     help="genomic region (FASTA with one sequence). Needed to run BLAST; optional with --blast")
    inp.add_argument("-b", "--blast", type=Path, metavar="FILE",
                     help="use an existing BLAST result instead of running BLAST. "
                          f"Must be -outfmt '{BLAST_OUTFMT}'")
    p.add_argument("-o", "--outdir", type=Path, default=Path("tast_results"), metavar="DIR",
                   help="output folder (default: tast_results)")

    flt = p.add_argument_group("filters")
    flt.add_argument("--strand", choices=STRAND_MODES, default="dominant",
                     help="alignment strand(s) kept per transcript. dominant: the strand with the highest "
                          "total bitscore; both: all alignments (TAST 1.0 behaviour) (default: dominant)")
    flt.add_argument("--max-intron", type=int, default=20000, metavar="BP",
                     help="remove transcripts whose largest gap between alignments is this long or longer "
                          "(default: 20000)")
    flt.add_argument("--min-coverage", type=float, default=25, metavar="PERCENT",
                     help="remove transcripts with this %% or less of their length aligned (default: 25)")
    flt.add_argument("--barcode-len", type=int, default=0, metavar="BP",
                     help="bases to trim from both ends of each transcript before finding ORFs, e.g. 39 for "
                          "the example data (default: 0)")
    flt.add_argument("--min-orf-len", type=int, default=26, metavar="NT",
                     help="ORFs must be longer than this; transcripts without one are not plotted (default: 26)")

    bl = p.add_argument_group("BLAST (only used when --blast is not given)")
    bl.add_argument("--evalue", type=float, default=1e-6, help="BLAST e-value cut-off (default: 1e-6)")
    bl.add_argument("--word-size", type=int, default=11, help="BLAST word size (default: 11)")
    bl.add_argument("--threads", type=int, default=1, help="CPU threads for BLAST (default: 1)")
    bl.add_argument("--blast-options", default="", metavar='"..."',
                    help='any other blastn options, in quotes, e.g. --blast-options "-dust no"')

    out = p.add_argument_group("plot")
    out.add_argument("--format", choices=["png", "pdf", "svg"], default="png", help="plot file type (default: png)")
    out.add_argument("--dpi", type=int, default=140, help="PNG resolution (default: 140)")
    out.add_argument("--numbering", choices=["transcript", "genome"], default="transcript",
                     help="how alignment regions are numbered on the plot. transcript: by position in the "
                          "transcript (0 = nearest the 5' end), so reverse-strand genes read ...2 1 0 and "
                          "rearrangements show as out-of-order numbers; genome: the order of TAST 1.0 "
                          "(default: transcript)")
    out.add_argument("--min-gap", type=int, default=0, metavar="BP",
                     help="minimum distance between transcripts sharing a plot row (default: 0)")
    out.add_argument("--no-row-numbers", action="store_true", help="do not number the plot rows")

    p.add_argument("-q", "--quiet", action="store_true", help="only print warnings and errors")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    started = time.time()

    def log(msg: str) -> None:
        if not args.quiet or msg.startswith("WARNING"):
            print(msg, file=sys.stderr, flush=True)

    try:
        return _run(args, log, started)
    except InputError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("Stopped.", file=sys.stderr)
        return 130


def _run(args, log, started) -> int:
    if args.blast is None and args.genome is None:
        raise InputError("give the genomic region with -g/--genome (to run BLAST), "
                         "or an existing BLAST result with -b/--blast")
    if not args.transcripts.is_file():
        raise InputError(f"transcripts FASTA not found: {args.transcripts}")
    if args.barcode_len < 0 or args.threads < 1:
        raise InputError("--barcode-len must be 0 or more and --threads 1 or more")
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)

    genome_id, seq_len = None, None
    if args.genome is not None:
        genome_id, genome_seq = read_genome(args.genome)
        seq_len = len(genome_seq)
        log(f"Genomic region {genome_id}: {seq_len:,} bp")

    #---- alignments
    blast_cmd = None
    if args.blast is None:
        from .blast import run_blast
        blast_path = out / "blast_hits.tsv"
        log("Aligning transcripts to the genomic region with BLAST (this can take a while)...")
        blast_cmd = run_blast(args.transcripts, args.genome, blast_path, evalue=args.evalue,
                              word_size=args.word_size, threads=args.threads,
                              extra_options=args.blast_options, log=log)
        log(f"BLAST result saved to {blast_path} (reuse it with --blast to skip BLAST next time)")
    else:
        blast_path = args.blast
    blast = read_blast(blast_path)
    if blast.empty:
        raise InputError("BLAST found no alignments between the transcripts and the genomic region.")
    subjects = blast["sseqid"].unique()
    if len(subjects) > 1:
        raise InputError(f"the BLAST result contains alignments to {len(subjects)} different sequences "
                         f"({', '.join(map(str, subjects[:3]))}...); TAST needs alignments to one genomic region")
    slen = int(blast["slen"].iloc[0])
    if seq_len is None:
        seq_len = slen
    elif seq_len != slen:
        log(f"WARNING: genome FASTA is {seq_len:,} bp but the BLAST subject is {slen:,} bp; "
            "check that BLAST was run against this sequence")
    log(f"{blast.shape[0]:,} alignments from {blast['qseqid'].nunique():,} transcripts")

    #---- analysis
    settings = Settings(max_intron=args.max_intron, min_coverage=args.min_coverage, strand=args.strand,
                        barcode_len=args.barcode_len, min_orf_len=args.min_orf_len, min_gap=args.min_gap,
                        numbering=args.numbering)
    result = run(blast, args.transcripts, out / "filtered_transcripts.fasta", settings, log=log)
    n_plotted = sum(len(r.transcripts) for r in result.rows)
    log(f"{n_plotted} transcripts arranged in {len(result.rows)} plot rows")

    #---- outputs
    write_lines(result.filtered_ids, out / "filtered_ids.txt")
    summary = result.summary.copy()
    summary["coverage"] = summary["coverage"].round(2)
    summary.to_csv(out / "transcripts.tsv", sep="\t")
    plotted = summary[summary["status"] == "plotted"].copy()
    plotted["genome_start"] = None
    plotted["genome_end"] = None
    for i, row in enumerate(result.rows):
        for name, (s, e) in zip(row.transcripts, row.ranges):
            plotted.loc[name, ["genome_start", "genome_end"]] = [s, e]
    plotted = plotted.reset_index().sort_values(["plot_row", "genome_start"])
    plotted[["plot_row", "qseqid", "genome_start", "genome_end", "n_alignments", "strand", "coverage",
             "max_intron", "largest_orf", "n_orfs", "region_order", "region_numbers"]].to_csv(out / "plot_rows.tsv", sep="\t", index=False)

    from .plot import draw, save
    plot_path = out / f"tast_plot.{args.format}"
    log("Drawing the plot...")
    fig = draw(result.rows, seq_len, row_numbers=not args.no_row_numbers, dpi=args.dpi)
    save(fig, plot_path, args.dpi, log=log)

    with open(out / "run_info.txt", "w", newline="\n") as f:
        f.write(f"lrRNAseq TAST {__version__}\n")
        f.write(f"date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("command: tast " + " ".join(sys.argv[1:]) + "\n")
        for k, v in sorted(vars(args).items()):
            f.write(f"{k}: {v}\n")
        if blast_cmd:
            f.write("blast command: " + " ".join(blast_cmd) + "\n")
        for w in result.warnings:
            f.write(f"warning: {w}\n")

    log(f"Done in {time.time() - started:.0f} s. Plot: {plot_path}  Row contents: {out / 'plot_rows.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

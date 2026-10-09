"""End-to-end tests on the example data in example/."""
import dataclasses
import json
import shutil
from pathlib import Path

import pytest

from tast.cli import main
from tast.files import read_blast
from tast.pipeline import Settings, run

ROOT = Path(__file__).resolve().parent.parent
EX = ROOT / "example"
BLAST = EX / "lrRNAseq_to_Ef318A10.out"
TRANSCRIPTS = EX / "all_transcripts.fasta"
GENOME = EX / "Ef318A10.fst"
EXPECTED = json.loads((Path(__file__).parent / "expected_example_rows.json").read_text())


@pytest.fixture(scope="module")
def blast():
    return read_blast(BLAST)


@pytest.mark.parametrize("strand", ["dominant", "both"])
@pytest.mark.parametrize("numbering", ["genome", "transcript"])
def test_rows_match_expected(blast, strand, numbering, tmp_path):
    #with numbering="genome" the rows are identical to the original notebook
    settings = Settings(strand=strand, barcode_len=39, numbering=numbering)
    res = run(blast, TRANSCRIPTS, tmp_path / "f.fasta", settings, log=lambda m: None)
    rows = json.loads(json.dumps([dataclasses.asdict(r) for r in res.rows], default=int))
    assert rows == EXPECTED[f"{strand}_{numbering}"]


def test_transcript_numbering_follows_strand(blast, tmp_path):
    #forward-strand transcripts read 0, 1, 2 ... along the genome; reverse-strand ones ... 2, 1, 0
    res = run(blast, TRANSCRIPTS, tmp_path / "f.fasta", Settings(barcode_len=39), log=lambda m: None)
    plotted = res.summary[res.summary["status"] == "plotted"]
    assert set(plotted.loc[plotted["strand"] == "+", "region_order"]) <= {"forward", "single", "mixed"}
    assert set(plotted.loc[plotted["strand"] == "-", "region_order"]) <= {"reverse", "single", "mixed"}
    assert (plotted["region_order"] == "forward").sum() == 45
    assert (plotted["region_order"] == "reverse").sum() == 4


def test_first_blast_line_is_kept(blast):
    #the BLAST file has no header line, so its first alignment must be read
    assert blast["qseqid"].iloc[0] == "m64045e_240201_155513/11602711/ccs_gut"
    assert len(blast) == 88482


def test_cli_with_blast_file(tmp_path):
    out = tmp_path / "res"
    code = main(["-t", str(TRANSCRIPTS), "-g", str(GENOME), "-b", str(BLAST), "-o", str(out),
                 "--barcode-len", "39", "-q"])
    assert code == 0
    for name in ["tast_plot.png", "plot_rows.tsv", "transcripts.tsv", "filtered_ids.txt",
                 "filtered_transcripts.fasta", "run_info.txt"]:
        assert (out / name).stat().st_size > 0
    assert len((out / "filtered_ids.txt").read_text().split()) == 56
    assert (out / "plot_rows.tsv").read_text().count("\n") == 57  #header + 56 transcripts


def test_cli_errors(tmp_path):
    assert main(["-t", str(TRANSCRIPTS), "-q"]) == 1                      #no genome or BLAST
    assert main(["-t", str(tmp_path / "missing.fa"), "-b", str(BLAST), "-q"]) == 1
    assert main(["-t", str(TRANSCRIPTS), "-b", str(BLAST), "-o", str(tmp_path / "x"),
                 "--min-coverage", "100", "-q"]) == 1                     #nothing passes


@pytest.mark.skipif(shutil.which("blastn") is None, reason="BLAST+ not installed")
@pytest.mark.parametrize("threads", [1, 2])
def test_cli_runs_blast(tmp_path, threads):
    out = tmp_path / f"blast{threads}"
    assert main(["-t", str(TRANSCRIPTS), "-g", str(GENOME), "-o", str(out), "--threads", str(threads), "-q"]) == 0
    #default BLAST settings reproduce the example BLAST file for these transcripts
    ids = {line[1:].split()[0] for line in TRANSCRIPTS.read_text().splitlines() if line.startswith(">")}
    original = {tuple(l.split("\t")[:8]) for l in BLAST.read_text().splitlines() if l.split("\t")[0] in ids}
    new = {tuple(l.split("\t")[:8]) for l in (out / "blast_hits.tsv").read_text().splitlines()}
    assert new == original


@pytest.mark.parametrize("extra", [[], ["--min-coverage", "85"]])  #15 rows, and a single row
def test_plot_fits_in_picture(tmp_path, extra):
    #nothing drawn (axis numbers, row numbers) may touch the edge of the saved picture
    import matplotlib.image as mpimg
    out = tmp_path / "res"
    assert main(["-t", str(TRANSCRIPTS), "-b", str(BLAST), "-o", str(out), "--barcode-len", "39", "-q"] + extra) == 0
    img = mpimg.imread(out / "tast_plot.png")[:, :, :3]
    for edge in (img[:5], img[-5:], img[:, :5], img[:, -5:]):
        assert edge.min() > 0.95

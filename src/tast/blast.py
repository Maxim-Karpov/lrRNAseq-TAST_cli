"""Running NCBI BLAST+ (blastn) to align transcripts to the genomic region."""
from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Callable, List, Optional, Tuple

from .files import BLAST_OUTFMT, InputError

BLAST_HELP = ("BLAST+ was not found. Install it from "
              "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ "
              "(on Windows use the 'win64.exe' installer) or with conda: "
              "conda install -c bioconda blast\n"
              "Alternatively, run BLAST yourself and pass the result with --blast.")


def _find(program: str) -> str:
    path = shutil.which(program)
    if path is None:
        raise InputError(BLAST_HELP)
    return path


def blast_version(blastn: str) -> Tuple[int, ...]:
    out = subprocess.run([blastn, "-version"], capture_output=True, text=True).stdout
    match = re.search(r"(\d+)\.(\d+)\.(\d+)", out)
    return tuple(int(x) for x in match.groups()) if match else (0, 0, 0)


def _run(cmd: List[str], cwd: Optional[str] = None) -> None:
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
    if proc.returncode != 0:
        raise InputError(f"BLAST failed (exit code {proc.returncode}):\n"
                         f"{' '.join(cmd)}\n{proc.stderr.strip()}")


def run_blast(transcripts: Path, genome: Path, out_path: Path, evalue: float = 1e-6,
              word_size: int = 11, task: str = "blastn", threads: int = 1,
              extra_options: str = "", log: Callable[[str], None] = print) -> List[str]:
    """Align transcripts (query) to the genomic region (subject) with blastn, writing the
    tabular output TAST needs to `out_path`. Returns the blastn command used.

    With one thread the genome is passed as -subject. With more threads a temporary BLAST
    database is built from the genome so that blastn can split the transcripts between
    threads (-mt_mode 1, BLAST+ 2.12 or newer).
    """
    blastn = _find("blastn")
    transcripts, genome, out_path = (Path(p).resolve() for p in (transcripts, genome, out_path))
    cmd = [blastn, "-task", task, "-query", str(transcripts), "-outfmt", BLAST_OUTFMT,
           "-evalue", str(evalue), "-word_size", str(word_size), "-out", str(out_path)]
    if extra_options:
        cmd += shlex.split(extra_options, posix=os.name != "nt")

    if threads <= 1:
        cmd += ["-subject", str(genome)]
        log("Running " + " ".join(Path(c).name if i == 0 else c for i, c in enumerate(cmd)))
        _run(cmd)
        return cmd

    makeblastdb = _find("makeblastdb")
    version = blast_version(blastn)
    with tempfile.TemporaryDirectory(prefix="tast_blastdb_") as tmp:
        #makeblastdb and BLAST database paths must not contain spaces, so the genome is copied
        #into the temporary folder and the database is referred to by a relative name
        shutil.copyfile(genome, Path(tmp) / "genome.fa")
        _run([makeblastdb, "-in", "genome.fa", "-dbtype", "nucl", "-out", "genome"], cwd=tmp)
        cmd += ["-db", "genome", "-num_threads", str(threads)]
        if version >= (2, 12, 0):
            cmd += ["-mt_mode", "1"]
        else:
            log(f"Note: BLAST+ {'.'.join(map(str, version))} cannot split queries between threads "
                "(needs 2.12 or newer), so --threads will have little effect")
        log("Running " + " ".join(Path(c).name if i == 0 else c for i, c in enumerate(cmd)))
        _run(cmd, cwd=tmp)
    return cmd

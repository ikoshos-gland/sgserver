"""Alignment: uses SnapGene's own bundled programs."""

from __future__ import annotations

import tempfile
from pathlib import Path

from Bio import SeqIO

from .core import SgError, run_tool

ALIGNERS = {
    "clustalo": "clustalOmega/clustalo",
    "muscle": "muscle/muscle",
}


def align(records, program: str = "clustalo", out: Path | None = None) -> str:
    """Run a multiple sequence alignment and return the aligned FASTA."""
    if program not in ALIGNERS:
        raise SgError(
            f"Unknown aligner: {program}. Options: {', '.join(ALIGNERS)}"
        )
    if len(records) < 2:
        raise SgError("At least 2 sequences are required for alignment.")

    with tempfile.TemporaryDirectory() as td:
        src = Path(td) / "in.fasta"
        dst = Path(td) / "out.fasta"
        SeqIO.write(records, str(src), "fasta")

        if program == "clustalo":
            args = ["-i", str(src), "-o", str(dst), "--outfmt=fasta", "--force"]
        else:  # muscle 5
            args = ["-align", str(src), "-output", str(dst)]

        r = run_tool(ALIGNERS[program], args, timeout=600)
        if not dst.exists():
            raise SgError(
                f"{program} could not align:\n{(r.stderr or r.stdout).strip()[:500]}"
            )
        text = dst.read_text()

    if out:
        out = Path(out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text)
    return text


def pairwise(a: str, b: str, mode: str = "global") -> str:
    """Align two sequences (parasail)."""
    with tempfile.TemporaryDirectory() as td:
        qf, df = Path(td) / "q.fa", Path(td) / "d.fa"
        qf.write_text(f">query\n{a}\n")
        df.write_text(f">subject\n{b}\n")
        flag = {"global": "nw", "local": "sw", "semi": "sg"}.get(mode, "nw")
        r = run_tool(
            "parasail/parasail_aligner",
            ["-a", f"{flag}_trace_striped_32", "-f", str(df), "-q", str(qf), "-O", "EMBOSS"],
            timeout=300,
        )
        if r.returncode != 0 and not r.stdout.strip():
            raise SgError(f"parasail failed:\n{r.stderr.strip()[:500]}")
        return r.stdout


def consensus_from_traces(fasta_paths: list[str]) -> str:
    """Assemble reads into a contig with CAP3."""
    with tempfile.TemporaryDirectory() as td:
        merged = Path(td) / "reads.fasta"
        recs = []
        for p in fasta_paths:
            recs.extend(SeqIO.parse(p, "fasta"))
        if not recs:
            raise SgError("No reads found.")
        SeqIO.write(recs, str(merged), "fasta")
        r = run_tool("cap3/cap3", [str(merged)], timeout=600)
        contigs = Path(str(merged) + ".cap.contigs")
        if not contigs.exists():
            raise SgError(f"CAP3 could not produce a contig:\n{r.stderr.strip()[:400]}")
        return contigs.read_text()

"""Align a sequencing read to a reference and report differences (SnapGene: Align to Reference)."""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner

from .core import SgError, revcomp

TRACE_FORMATS = {".ab1": "abi", ".abi": "abi"}


def read_trace(path: str | Path) -> str:
    """Read the base sequence from a sequencing file (.ab1). Plain sequence files work too."""
    path = Path(path)
    if not path.exists():
        raise SgError(f"File not found: {path}")
    suffix = path.suffix.lower()
    if suffix in TRACE_FORMATS:
        rec = SeqIO.read(str(path), TRACE_FORMATS[suffix])
        return str(rec.seq).upper()
    if suffix == ".scf":
        raise SgError(
            "SCF trace files cannot be read (Biopython does not support them).\n"
            "Open the file in SnapGene, export it as FASTA, then pass that instead."
        )
    from .core import load_many

    return str(load_many(path)[0].seq).upper()


def _aligner() -> PairwiseAligner:
    a = PairwiseAligner()
    a.mode = "local"
    a.match_score = 2
    a.mismatch_score = -3
    a.open_gap_score = -7
    a.extend_gap_score = -2
    return a


def align_read(reference, read: str, circular_aware: bool = True) -> dict:
    """Align the read to the reference; return best strand, position and differences."""
    read = read.upper().strip()
    if len(read) < 20:
        raise SgError(f"Read too short: {len(read)} bp")
    ref = str(reference.seq).upper()
    n = len(ref)
    if circular_aware and getattr(reference, "circular", False):
        ref_search = ref + ref[: min(len(read), n)]
    else:
        ref_search = ref

    aligner = _aligner()
    best = None
    for strand, probe in ((1, read), (-1, revcomp(read))):
        alignment = aligner.align(ref_search, probe)[0]
        if best is None or alignment.score > best[0].score:
            best = (alignment, strand, probe)
    alignment, strand, probe = best

    ref_aln, read_aln = alignment[0], alignment[1]
    ref_start = int(alignment.coordinates[0][0])

    differences = []
    ri = ref_start
    for r_char, q_char in zip(ref_aln, read_aln):
        if r_char == "-":
            differences.append({"position": (ri % n) + 1, "type": "insertion",
                                "reference": "-", "read": q_char})
            continue
        if q_char == "-":
            differences.append({"position": (ri % n) + 1, "type": "deletion",
                                "reference": r_char, "read": "-"})
        elif r_char != q_char:
            differences.append({"position": (ri % n) + 1, "type": "substitution",
                                "reference": r_char, "read": q_char})
        ri += 1

    aligned_len = sum(1 for a, b in zip(ref_aln, read_aln) if a != "-" and b != "-")
    matches = sum(1 for a, b in zip(ref_aln, read_aln) if a == b and a != "-")
    return {
        "strand": strand,
        "start": (ref_start % n) + 1,
        "end": ((ref_start + aligned_len - 1) % n) + 1,
        "aligned": aligned_len,
        "read_length": len(read),
        "identity_pct": 100.0 * matches / aligned_len if aligned_len else 0.0,
        "differences": differences,
        "score": alignment.score,
    }


def format_result(res: dict, max_diff: int = 25) -> str:
    direction = "forward (+)" if res["strand"] == 1 else "reverse (−)"
    lines = [
        f"Read of {res['read_length']} bp aligned on the {direction} strand",
        f"Reference {res['start']}..{res['end']} "
        f"({res['aligned']} bp overlap, identity {res['identity_pct']:.2f}%)",
    ]
    if not res["differences"]:
        lines.append("No differences — the read matches the reference exactly.")
        return "\n".join(lines)
    lines.append(f"\n{len(res['differences'])} differences:")
    for d in res["differences"][:max_diff]:
        lines.append(f"  {d['position']:>6}  {d['type']:<12} {d['reference']} → {d['read']}")
    if len(res["differences"]) > max_diff:
        lines.append(f"  ... and {len(res['differences']) - max_diff} more")
    return "\n".join(lines)

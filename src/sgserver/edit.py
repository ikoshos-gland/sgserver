"""Sequence editing: mutagenesis, origin, flip, feature annotation."""

from __future__ import annotations

import re

from Bio.SeqFeature import SeqFeature, SimpleLocation
from pydna.dseqrecord import Dseqrecord

from .core import SgError, revcomp

VALID = set("ACGTRYSWKMBDHVN")


def _check(seq: str, what: str) -> str:
    s = seq.upper().replace(" ", "")
    bad = set(s) - VALID
    if bad:
        raise SgError(f"Invalid letter(s) in {what}: {', '.join(sorted(bad))}")
    return s


def _rebuild(rec, new_seq: str, change_at: int, delta: int):
    """Build a record from the new sequence; keep/shift features where possible."""
    out = Dseqrecord(new_seq, circular=bool(getattr(rec, "circular", False)))
    out.name = rec.name
    out.id = rec.id
    out.description = rec.description
    for f in rec.features:
        start, end = int(f.location.start), int(f.location.end)
        if delta == 0 or end <= change_at:
            # if the length does not change (point mutation) no feature shifts
            ns, ne = start, end
        elif start >= change_at - min(delta, 0):
            ns, ne = start + delta, end + delta
        else:
            continue  # the change passes through the feature: moving it would be misleading
        if ns < 0 or ne > len(new_seq) or ns >= ne:
            continue
        out.features.append(
            SeqFeature(
                location=SimpleLocation(ns, ne, strand=f.location.strand),
                type=f.type,
                qualifiers=dict(f.qualifiers),
            )
        )
    return out


def substitute(rec, start: int, new: str):
    """Replace with `new` starting at 1-based position `start` (same length)."""
    new = _check(new, "new sequence")
    seq = str(rec.seq)
    n = len(seq)
    if not (1 <= start <= n):
        raise SgError(f"Position out of range: {start} (sequence is {n} bp)")
    end = start - 1 + len(new)
    if end > n and not getattr(rec, "circular", False):
        raise SgError(f"Change runs past the end of the sequence: {start}..{end} > {n}")
    if end <= n:
        out_seq = seq[: start - 1] + new + seq[end:]
    else:  # wraps around the circle
        over = end - n
        out_seq = new[len(new) - over :] + seq[over : start - 1] + new[: len(new) - over]
    return _rebuild(rec, out_seq, start - 1, 0)


def insert(rec, after: int, fragment: str):
    """Insert `fragment` AFTER 1-based position `after`. after=0 → at the very start."""
    fragment = _check(fragment, "inserted sequence")
    seq = str(rec.seq)
    n = len(seq)
    if not (0 <= after <= n):
        raise SgError(f"Position out of range: {after} (sequence is {n} bp)")
    out_seq = seq[:after] + fragment + seq[after:]
    return _rebuild(rec, out_seq, after, len(fragment))


def delete(rec, start: int, end: int):
    """Delete the 1-based, inclusive range [start, end]."""
    seq = str(rec.seq)
    n = len(seq)
    if not (1 <= start <= end <= n):
        raise SgError(f"Invalid range: {start}..{end} (sequence is {n} bp)")
    out_seq = seq[: start - 1] + seq[end:]
    if not out_seq:
        raise SgError("Cannot delete the entire sequence.")
    return _rebuild(rec, out_seq, start - 1, -(end - start + 1))


def set_origin(rec, position: int):
    """Move the origin of a circular plasmid to 1-based `position`."""
    if not getattr(rec, "circular", False):
        raise SgError("The origin can only be changed on circular sequences.")
    n = len(rec.seq)
    if not (1 <= position <= n):
        raise SgError(f"Position out of range: {position} (sequence is {n} bp)")
    out = rec.shifted(position - 1)
    out.name, out.id = rec.name, rec.id
    return out


def flip(rec):
    """Take the reverse complement (SnapGene: 'Flip Sequence')."""
    out = rec.reverse_complement()
    out.name, out.id = rec.name, rec.id
    out.description = rec.description
    return out


def add_feature(rec, start: int, end: int, label: str, kind: str = "misc_feature",
                strand: int = 1):
    """Add a feature over the 1-based, inclusive range."""
    n = len(rec.seq)
    if not (1 <= start <= end <= n):
        raise SgError(f"Invalid range: {start}..{end} (sequence is {n} bp)")
    rec.features.append(
        SeqFeature(
            location=SimpleLocation(start - 1, end, strand=strand),
            type=kind,
            qualifiers={"label": [label]},
        )
    )
    rec.features.sort(key=lambda f: int(f.location.start))
    return rec


# --------------------------------------------------------------------------
# site-directed mutagenesis primers
# --------------------------------------------------------------------------


def mutagenesis_primers(rec, position: int, new: str, flank: int = 15):
    """QuikChange-style mutagenesis primer pair, complementary to each other."""
    new = _check(new, "new sequence")
    seq = str(rec.seq).upper()
    n = len(seq)
    if not (1 <= position <= n):
        raise SgError(f"Position out of range: {position} (sequence is {n} bp)")
    circular = bool(getattr(rec, "circular", False))
    doubled = seq + seq if circular else seq

    left_start = position - 1 - flank
    if left_start < 0:
        if not circular:
            raise SgError(f"Too close to the edge: must be at least {flank} bp from the start.")
        left_start += n
    left = doubled[left_start : left_start + flank]
    right_start = (position - 1 + len(new)) % n if circular else position - 1 + len(new)
    right = doubled[right_start : right_start + flank]
    if len(left) < flank or len(right) < flank:
        raise SgError("Not enough flanking sequence; reduce `flank`.")

    fwd = left + new + right
    return {"fwd": fwd, "rev": revcomp(fwd), "old": seq[position - 1 : position - 1 + len(new)],
            "new": new, "position": position}


# --------------------------------------------------------------------------
# CRISPR rehber dizileri
# --------------------------------------------------------------------------


def find_guides(rec, pam: str = "NGG", length: int = 20, limit: int | None = None):
    """SpCas9 guide candidates: protospacer + PAM, on both strands.

    Returns: [{start, end, strand, guide, pam, gc}] (1-based, inclusive)
    """
    pam = pam.upper()
    iupac = {"A": "A", "C": "C", "G": "G", "T": "T", "N": "[ACGT]",
             "R": "[AG]", "Y": "[CT]", "W": "[AT]", "S": "[GC]"}
    try:
        pam_rx = "".join(iupac[c] for c in pam)
    except KeyError as e:
        raise SgError(f"Invalid letter in PAM: {e.args[0]}") from e

    seq = str(rec.seq).upper()
    n = len(seq)
    circular = bool(getattr(rec, "circular", False))
    pad = length + len(pam) if circular else 0
    hay = seq + seq[:pad]

    out = []
    for strand in (1, -1):
        target = hay if strand == 1 else revcomp(hay)
        for m in re.finditer(f"(?=({pam_rx}))", target):
            p = m.start()
            if p < length:
                continue
            guide = target[p - length : p]
            if len(guide) < length or "N" in guide:
                continue
            if strand == 1:
                start0 = p - length
            else:
                start0 = len(target) - p
            if start0 >= n:
                continue
            gc = 100.0 * (guide.count("G") + guide.count("C")) / length
            out.append({
                "start": start0 + 1,
                "end": start0 + length,
                "strand": strand,
                "guide": guide,
                "pam": target[p : p + len(pam)],
                "gc": gc,
            })
    out.sort(key=lambda g: (g["start"], -g["strand"]))
    return out[:limit] if limit else out

"""Sequence inspection: summary, restriction enzymes, virtual gel, ORFs/translation, search."""

from __future__ import annotations

import re
from pathlib import Path

from Bio.Restriction import AllEnzymes, CommOnly, RestrictionBatch
from Bio.Seq import Seq

from .core import SgError, fmt_bp, label_of, revcomp, topology_str

# --------------------------------------------------------------------------
# summary
# --------------------------------------------------------------------------


def gc_percent(seq: str) -> float:
    s = seq.upper()
    return 100.0 * (s.count("G") + s.count("C")) / len(s) if s else 0.0


def summarize(rec, show_features: int = 40) -> str:
    seq = str(rec.seq)
    lines = [
        f"Name      : {rec.name}",
        f"Length    : {fmt_bp(len(seq))}",
        f"Topology  : {topology_str(rec)}",
        f"GC content: {gc_percent(seq):.1f}%",
    ]
    desc = (rec.description or "").strip()
    if desc and desc != rec.name:
        lines.append(f"Description: {desc}")

    feats = [f for f in rec.features if f.type != "source"]
    lines.append(f"\nFeatures ({len(feats)}):")
    if not feats:
        lines.append("  (none)")
    for f in feats[:show_features]:
        start = int(f.location.start) + 1
        end = int(f.location.end)
        strand = {1: "+", -1: "-"}.get(f.location.strand, " ")
        lines.append(
            f"  {start:>6}..{end:<6} {strand}  {f.type:<14} {label_of(f)}"
        )
    if len(feats) > show_features:
        lines.append(f"  ... and {len(feats) - show_features} more")
    return "\n".join(lines)


# --------------------------------------------------------------------------
# restriction enzymes
# --------------------------------------------------------------------------


def _batch(names: list[str] | None, common_only: bool = True) -> RestrictionBatch:
    if not names:
        return RestrictionBatch(list(CommOnly)) if common_only else RestrictionBatch(list(AllEnzymes))
    out = []
    for n in names:
        enz = next((e for e in AllEnzymes if str(e) == n), None)
        if enz is None:
            raise SgError(f"Unknown enzyme: {n}")
        out.append(enz)
    return RestrictionBatch(out)


def analyze_sites(rec, enzymes: list[str] | None = None, max_cuts: int | None = None,
                  common_only: bool = True) -> list[tuple[str, list[int]]]:
    """Return the cut positions for each enzyme (1-based)."""
    batch = _batch(enzymes, common_only)
    circular = bool(getattr(rec, "circular", False))
    result = batch.search(Seq(str(rec.seq)), linear=not circular)
    rows = []
    for enz, sites in result.items():
        if not sites:
            continue
        if max_cuts is not None and len(sites) > max_cuts:
            continue
        rows.append((str(enz), sorted(sites)))
    rows.sort(key=lambda r: (len(r[1]), r[0]))
    return rows


def format_sites(rows, title: str = "Restriction sites") -> str:
    if not rows:
        return f"{title}: (none)"
    lines = [f"{title} ({len(rows)} enzymes):"]
    for name, sites in rows:
        pos = ", ".join(str(s) for s in sites[:12])
        more = f" ... (+{len(sites) - 12})" if len(sites) > 12 else ""
        lines.append(f"  {name:<12} {len(sites):>3} cuts   {pos}{more}")
    return "\n".join(lines)


def digest(rec, enzymes: list[str]):
    """Cut with the given enzymes and return the fragments."""
    batch = _batch(enzymes, common_only=False)
    frags = rec.cut(*list(batch))
    if not frags:
        raise SgError(
            f"{', '.join(enzymes)} does not cut this sequence (no sites)."
        )
    return list(frags)


def fragment_size(frag) -> int:
    """Standard size (bp) of a digest fragment.

    pydna's `len()` counts sticky ends twice: a 7235 bp circle cut with BamHI
    comes out as 7239, and the fragments do not add up to the plasmid size.
    The size used in restriction maps is the distance between the cut points
    on the top strand; these always sum to the total length.
    """
    try:
        return len(frag.seq.watson)
    except AttributeError:
        return len(frag)


def _end_label(kind: str, overhang: str) -> str:
    if kind == "blunt":
        return "blunt"
    side = "5' overhang" if kind.startswith("5") else "3' overhang"
    return f"{side} {overhang.upper()}" if overhang else side


def format_digest(frags, enzymes: list[str]) -> str:
    lines = [f"{' + '.join(enzymes)} digest → {len(frags)} fragments:"]
    for i, f in enumerate(sorted(frags, key=fragment_size, reverse=True), 1):
        ends = ""
        try:
            kind5, ov5 = f.seq.five_prime_end()
            kind3, ov3 = f.seq.three_prime_end()
            ends = f"   ends: {_end_label(kind5, ov5)} / {_end_label(kind3, ov3)}"
        except Exception:
            pass
        lines.append(f"  {i:>2}. {fmt_bp(fragment_size(f)):>12}{ends}")
    return "\n".join(lines)


# --------------------------------------------------------------------------
# virtual gel
# --------------------------------------------------------------------------

LADDER_1KB = [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 500, 250]


def gel_image(lanes: dict[str, list[int]], out: Path, ladder: list[int] | None = None) -> Path:
    """Draw a simple virtual agarose gel. lanes: {'lane name': [band sizes]}"""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    ladder = ladder or LADDER_1KB
    all_lanes = {"Ladder": ladder, **lanes}
    if not any(s > 0 for name, v in all_lanes.items() if name != "Ladder" for s in v):
        raise SgError("No bands to show on the gel (all lanes except the ladder are empty).")
    sizes = [s for v in all_lanes.values() for s in v if s > 0]
    top, bottom = max(sizes) * 1.6, max(min(sizes) * 0.6, 20)

    def y(size):  # agarose: migration is inversely related to log(size)
        return (np.log10(top) - np.log10(size)) / (np.log10(top) - np.log10(bottom))

    n = len(all_lanes)
    fig, ax = plt.subplots(figsize=(1.4 * n + 1.5, 7))
    ax.set_facecolor("#1a1a1a")
    fig.patch.set_facecolor("#1a1a1a")

    for i, (name, bands) in enumerate(all_lanes.items()):
        for b in bands:
            if b <= 0:
                continue
            yy = y(b)
            ax.add_patch(
                plt.Rectangle((i + 0.12, yy - 0.006), 0.76, 0.012,
                              color="#e8e8e8", alpha=0.95)
            )
            if name == "Ladder":
                ax.text(i + 0.05, yy, f"{b}", ha="right", va="center",
                        color="#cccccc", fontsize=7)
            else:
                ax.text(i + 0.94, yy, f"{b}", ha="left", va="center",
                        color="#8fd98f", fontsize=7)
        ax.text(i + 0.5, -0.04, name, ha="center", va="bottom",
                color="#ffffff", fontsize=9, transform=ax.get_xaxis_transform())

    ax.set_xlim(-0.6, n + 0.4)
    ax.set_ylim(1.03, -0.08)
    ax.axis("off")
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return out


# --------------------------------------------------------------------------
# ORF / translation
# --------------------------------------------------------------------------


def find_orfs(rec, min_aa: int = 100, table: int = 1):
    """Find ORFs in 3 frames on both strands."""
    seq = str(rec.seq).upper()
    circular = bool(getattr(rec, "circular", False))
    search = seq + seq[: min(len(seq), 3 * min_aa * 3)] if circular else seq
    found = []
    for strand, s in ((1, search), (-1, revcomp(search))):
        for frame in range(3):
            i = frame
            while i < len(s) - 2:
                if s[i : i + 3] == "ATG":
                    j = i
                    while j < len(s) - 2:
                        codon = s[j : j + 3]
                        if codon in ("TAA", "TAG", "TGA"):
                            aa_len = (j - i) // 3
                            if aa_len >= min_aa:
                                if strand == 1:
                                    start, end = i, j + 3
                                else:
                                    start, end = len(s) - (j + 3), len(s) - i
                                if start < len(seq):
                                    found.append((start + 1, min(end, len(seq) * 2),
                                                  strand, aa_len))
                            i = j
                            break
                        j += 3
                i += 3
    # keep the longest of nested ORFs that share the same end
    best: dict = {}
    for start, end, strand, aa in found:
        key = (end, strand)
        if key not in best or aa > best[key][3]:
            best[key] = (start, end, strand, aa)
    return sorted(best.values(), key=lambda r: -r[3])


def translate_region(rec, start: int, end: int, strand: int = 1, table: int = 1) -> str:
    """Translate a 1-based, inclusive range to protein."""
    seq = str(rec.seq)
    if getattr(rec, "circular", False) and end > len(seq):
        seq = seq + seq
    sub = seq[start - 1 : end]
    if strand == -1:
        sub = revcomp(sub)
    sub = sub[: len(sub) - len(sub) % 3]
    return str(Seq(sub).translate(table=table))


# --------------------------------------------------------------------------
# search
# --------------------------------------------------------------------------

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": ".",
}


def find_sequence(rec, query: str, both_strands: bool = True):
    """Search for a sequence (IUPAC ambiguity codes supported). Returns a list of (start, end, strand)."""
    q = query.upper().replace(" ", "")
    bad = set(q) - set(IUPAC)
    if bad:
        raise SgError(f"Invalid letter(s): {', '.join(sorted(bad))}")
    seq = str(rec.seq).upper()
    if getattr(rec, "circular", False):
        seq = seq + seq[: len(q) - 1]
    targets = [(1, q)]
    if both_strands and revcomp(q) != q:
        targets.append((-1, revcomp(q)))
    hits = []
    for strand, pattern in targets:
        rx = re.compile("".join(IUPAC[c] for c in pattern))
        for m in rx.finditer(seq):
            pos = m.start() + 1
            if pos <= len(rec.seq):
                hits.append((pos, pos + len(q) - 1, strand))
    return sorted(set(hits))

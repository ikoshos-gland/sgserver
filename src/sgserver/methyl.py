"""Does Dam / Dcm methylation block restriction sites?

Plasmids propagated in E. coli carry two methylases:

* **Dam**  — methylates the A in GATC.
* **Dcm**  — methylates the second C in CCWGG (CCAGG or CCTGG).

If an enzyme's recognition site covers a methylated base (it overlaps GATC or
CCWGG), the enzyme may be unable to cut there. SnapGene maps mark such sites
with an asterisk (e.g. `ClaI*`, `XbaI*`).

What is computed here is "overlap": does the site include a methylated base?
Whether the enzyme is really sensitive to that methylation depends on the
enzyme (REBASE data), so the result is a "may be blocked" warning.
"""

from __future__ import annotations

import re

from Bio.Seq import Seq

DAM_SITE = "GATC"
DAM_METHYLATED_OFFSET = 1          # the A inside GATC
DCM_PATTERN = re.compile("CC[AT]GG")
DCM_METHYLATED_OFFSET = 1          # the second C inside CCWGG


def _methylated_positions(seq: str, circular: bool) -> dict[int, set[str]]:
    """Methylated base positions (0-based) → {'Dam', 'Dcm'}."""
    n = len(seq)
    pad = 5 if circular else 0
    hay = seq + seq[:pad] if circular else seq
    out: dict[int, set[str]] = {}

    def mark(pos: int, kind: str) -> None:
        out.setdefault(pos % n, set()).add(kind)

    # Dam: both strands are methylated; GATC is palindromic so one scan is enough
    start = hay.find(DAM_SITE)
    while start != -1:
        if start < n:
            mark(start + DAM_METHYLATED_OFFSET, "Dam")
            mark(start + 2, "Dam")  # the A on the opposite strand (GATC → GATC)
        start = hay.find(DAM_SITE, start + 1)

    # Dcm: CCWGG and its opposite strand
    for strand_seq, flip in ((hay, False), (str(Seq(hay).reverse_complement()), True)):
        for m in DCM_PATTERN.finditer(strand_seq):
            pos = m.start() + DCM_METHYLATED_OFFSET
            if flip:
                pos = len(strand_seq) - 1 - pos
            if 0 <= pos < n + pad:
                mark(pos, "Dcm")
    return out


# Well-established behaviours. Only textbook-level, verifiable cases go here;
# any doubtful enzyme is reported as "overlap".
KNOWN_DAM_BLOCKED = {
    "MboI",    # GATC — Dam methylation blocks cutting
    "DpnII",   # GATC — Dam methylation blocks cutting
}
KNOWN_DAM_REQUIRED = {
    "DpnI",    # cuts ONLY Dam-methylated GATC
}
# GGATCC contains GATC but BamHI is not affected by Dam: proof that overlap
# alone is not enough. Hence "overlap" and "blocked" are kept apart below.
KNOWN_DAM_INSENSITIVE = {"BamHI", "BglII", "Sau3AI", "PvuI", "KpnI"}


def classify(enzyme_name: str, kinds: set[str]) -> str:
    """Return 'blocked' | 'required' | 'unaffected' | 'overlap'."""
    if "Dam" in kinds:
        if enzyme_name in KNOWN_DAM_REQUIRED:
            return "required"
        if enzyme_name in KNOWN_DAM_BLOCKED:
            return "blocked"
        if enzyme_name in KNOWN_DAM_INSENSITIVE:
            return "unaffected"
    return "overlap"


def blocked_sites(rec, enzyme_sites: dict[str, list[int]], enzyme_lengths: dict[str, int]):
    """Which cut sites overlap methylation?

    enzyme_sites: {enzyme name: [1-based cut positions]}
    enzyme_lengths: {enzyme name: recognition site length}

    Returns: {enzyme name: {position: {'Dam'|'Dcm'}}}
    """
    seq = str(rec.seq).upper()
    n = len(seq)
    circular = bool(getattr(rec, "circular", False))
    methylated = _methylated_positions(seq, circular)
    if not methylated:
        return {}

    result: dict[str, dict[int, set[str]]] = {}
    for name, positions in enzyme_sites.items():
        site_len = enzyme_lengths.get(name, 6)
        for cut in positions:
            # Possible range of the recognition site around the cut point
            # (Bio.Restriction gives the cut position, not the site start.)
            lo = cut - site_len - 1
            hi = cut + site_len
            kinds: set[str] = set()
            for i in range(lo, hi):
                kinds |= methylated.get(i % n, set())
            if kinds:
                result.setdefault(name, {})[cut] = kinds
    return result


ETIKET = {
    "blocked": "BLOCKED",
    "required": "methylation REQUIRED",
    "unaffected": "unaffected",
    "overlap": "overlaps — check the enzyme's sensitivity",
}


def describe_blocked(blocked: dict) -> str:
    if not blocked:
        return "No cut sites overlap Dam/Dcm sites."
    lines = ["Cut sites overlapping Dam/Dcm methylation:"]
    for name in sorted(blocked):
        for cut, kinds in sorted(blocked[name].items()):
            durum = ETIKET[classify(name, kinds)]
            lines.append(f"  {name:<12} {cut:>6}   {'+'.join(sorted(kinds)):<8} {durum}")
    lines.append(
        "\nOverlap alone does not mean blocked: GGATCC contains GATC "
        "but BamHI is not affected by Dam.\nCheck the enzyme's sensitivity in REBASE or "
        "the supplier catalogue; if it is sensitive, use a dam-/dcm- strain."
    )
    return "\n".join(lines)

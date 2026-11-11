"""PCR and primer tasks: design, amplification, Tm, dimer/hairpin checks."""

from __future__ import annotations

import primer3
from pydna.amplify import pcr as _pcr
from pydna.primer import Primer

from .core import SgError, fmt_bp, revcomp, run_tool

# --------------------------------------------------------------------------
# Tm and primer quality
# --------------------------------------------------------------------------


def tm(seq: str, dna_nm: float = 500.0, mv: float = 50.0, dv: float = 1.5,
       dntp: float = 0.6) -> float:
    """Tm (°C) using SantaLucia thermodynamics. SnapGene uses the same approach."""
    return primer3.calc_tm(
        seq.upper(), mv_conc=mv, dv_conc=dv, dntp_conc=dntp, dna_conc=dna_nm
    )


def primer_report(seq: str) -> dict:
    s = seq.upper().strip()
    if not s:
        raise SgError("Empty primer.")
    hp = primer3.calc_hairpin(s)
    hd = primer3.calc_homodimer(s)
    gc = 100.0 * (s.count("G") + s.count("C")) / len(s)
    return {
        "sequence": s,
        "length": len(s),
        "tm": tm(s),
        "gc": gc,
        "hairpin_tm": hp.tm if hp.structure_found else None,
        "hairpin_dg": hp.dg / 1000 if hp.structure_found else None,
        "homodimer_tm": hd.tm if hd.structure_found else None,
        "homodimer_dg": hd.dg / 1000 if hd.structure_found else None,
        "gc_clamp_3p": s[-5:].count("G") + s[-5:].count("C"),
    }


def heterodimer(a: str, b: str) -> dict:
    r = primer3.calc_heterodimer(a.upper(), b.upper())
    return {
        "tm": r.tm if r.structure_found else None,
        "dg": r.dg / 1000 if r.structure_found else None,
    }


def format_primer_report(rep: dict, name: str = "Primer") -> str:
    lines = [
        f"{name}",
        f"  Sequence    : {rep['sequence']}",
        f"  Length      : {rep['length']} nt",
        f"  Tm          : {rep['tm']:.1f} °C",
        f"  GC          : %{rep['gc']:.1f}",
        f"  3' GC clamp : {rep['gc_clamp_3p']}/5",
    ]
    if rep["hairpin_tm"] is not None:
        warn = "  ⚠ high" if rep["hairpin_tm"] > 45 else ""
        lines.append(f"  Hairpin     : Tm {rep['hairpin_tm']:.1f} °C, "
                     f"ΔG {rep['hairpin_dg']:.1f} kcal/mol{warn}")
    else:
        lines.append("  Hairpin     : none")
    if rep["homodimer_tm"] is not None:
        warn = "  ⚠ high" if rep["homodimer_tm"] > 45 else ""
        lines.append(f"  Self-dimer  : Tm {rep['homodimer_tm']:.1f} °C, "
                     f"ΔG {rep['homodimer_dg']:.1f} kcal/mol{warn}")
    else:
        lines.append("  Self-dimer  : none")
    return "\n".join(lines)


# --------------------------------------------------------------------------
# primer design
# --------------------------------------------------------------------------


def design_primers(
    rec,
    start: int,
    end: int,
    target_tm: float = 58.0,
    min_len: int = 18,
    max_len: int = 30,
    fwd_tail: str = "",
    rev_tail: str = "",
) -> dict:
    """Design a primer pair that amplifies the 1-based [start, end] region."""
    seq = str(rec.seq).upper()
    n = len(seq)
    if getattr(rec, "circular", False):
        seq = seq + seq
    if not (1 <= start <= n) or end < start:
        raise SgError(f"Invalid range: {start}..{end} (sequence is {n} bp)")

    region = seq[start - 1 : end]
    if len(region) < min_len * 2:
        raise SgError(f"Region too short: {len(region)} bp")

    def pick(template: str) -> str:
        best, best_d = None, 1e9
        for L in range(min_len, max_len + 1):
            cand = template[:L]
            d = abs(tm(cand) - target_tm)
            if d < best_d:
                best, best_d = cand, d
        return best

    fwd = pick(region)
    rev = pick(revcomp(region))

    return {
        "fwd": fwd_tail.upper() + fwd,
        "rev": rev_tail.upper() + rev,
        "fwd_anneal": fwd,
        "rev_anneal": rev,
        "product_bp": len(region),
        "region": (start, end),
    }


# --------------------------------------------------------------------------
# PCR
# --------------------------------------------------------------------------


def amplify(rec, fwd: str, rev: str, limit: int = 14,
            fwd_name: str = "fwd_primer", rev_name: str = "rev_primer"):
    """Compute the PCR product for a primer pair."""
    try:
        f = Primer(fwd.upper(), id=fwd_name, name=fwd_name)
        r = Primer(rev.upper(), id=rev_name, name=rev_name)
        product = _pcr(f, r, rec, limit=limit)
    except Exception as e:
        raise SgError(
            f"No PCR product found: {e}\n"
            "The primers may not anneal to the template. "
            "Try lowering the match length with --limit."
        ) from e
    return product


def format_amplicon(amp) -> str:
    lines = [
        f"PCR product: {fmt_bp(len(amp))}",
        f"  Forward primer: {amp.forward_primer.seq}  (Tm {tm(str(amp.forward_primer.seq)):.1f} °C)",
        f"  Reverse primer: {amp.reverse_primer.seq}  (Tm {tm(str(amp.reverse_primer.seq)):.1f} °C)",
    ]
    het = heterodimer(str(amp.forward_primer.seq), str(amp.reverse_primer.seq))
    if het["tm"] is not None:
        warn = "  ⚠" if het["tm"] > 45 else ""
        lines.append(f"  Primer dimer: Tm {het['tm']:.1f} °C, ΔG {het['dg']:.1f}{warn}")
    else:
        lines.append("  Primer dimer: none")
    return "\n".join(lines)


# --------------------------------------------------------------------------
# RNA/DNA folding (ViennaRNA bundled with SnapGene)
# --------------------------------------------------------------------------


def fold(seq: str) -> tuple[str, float]:
    """Secondary structure and free energy via RNAfold."""
    import subprocess

    from .core import TOOLS, tool_env

    p = subprocess.run(
        [str(TOOLS / "viennaRna" / "RNAfold"), "--noPS"],
        input=seq.upper() + "\n",
        capture_output=True,
        text=True,
        env=tool_env(),
        timeout=120,
    )
    out = p.stdout.strip().splitlines()
    if len(out) < 2:
        raise SgError(f"RNAfold failed: {p.stderr.strip()}")
    structure_line = out[-1]
    struct = structure_line.split(" ")[0]
    try:
        energy = float(structure_line.split("(")[-1].rstrip(")").strip())
    except ValueError:
        energy = float("nan")
    return struct, energy

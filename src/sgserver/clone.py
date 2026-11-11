"""Cloning: ligation, Gibson, Golden Gate, Gateway, and SnapGene history."""

from __future__ import annotations

from pathlib import Path

from Bio.Restriction import AllEnzymes
from pydna import assembly2 as a2
from pydna.dseqrecord import Dseqrecord

from .core import SgError, fmt_bp, load

METHODS = {
    "ligation": "ligation",
    "restriction": "restriction",
    "gibson": "gibson",
    "goldengate": "goldengate",
    "golden-gate": "goldengate",
    "infusion": "infusion",
    "in-fusion": "infusion",
    "gateway-bp": "gateway-bp",
    "gateway-lr": "gateway-lr",
}


def _enzymes(names: list[str]):
    out = []
    for n in names:
        e = next((x for x in AllEnzymes if str(x) == n), None)
        if e is None:
            raise SgError(f"Unknown enzyme: {n}")
        out.append(e)
    return out


def assemble(
    frags: list[Dseqrecord],
    method: str,
    enzymes: list[str] | None = None,
    limit: int = 25,
    circular_only: bool = True,
) -> list[Dseqrecord]:
    """Assemble the fragments with the given method and return the possible products."""
    m = METHODS.get(method.lower())
    if m is None:
        raise SgError(
            f"Unknown method: {method}\n"
            f"Options: {', '.join(sorted(set(METHODS.values())))}"
        )

    if m in ("restriction", "goldengate") and not enzymes:
        raise SgError(f"'{method}' requires at least one enzyme via --enzyme.")

    if m == "ligation":
        products = a2.ligation_assembly(frags, allow_blunt=True, circular_only=circular_only)
    elif m == "restriction":
        products = a2.restriction_ligation_assembly(
            frags, _enzymes(enzymes), circular_only=circular_only
        )
    elif m == "goldengate":
        products = a2.golden_gate_assembly(
            frags, _enzymes(enzymes), circular_only=circular_only
        )
    elif m == "gibson":
        products = a2.gibson_assembly(frags, limit=limit, circular_only=circular_only)
    elif m == "infusion":
        products = a2.in_fusion_assembly(frags, limit=limit, circular_only=circular_only)
    elif m.startswith("gateway"):
        products = a2.gateway_assembly(
            frags, reaction_type=m.split("-")[1].upper(), circular_only=circular_only
        )
    else:  # pragma: no cover
        raise SgError(f"Unsupported method: {m}")

    if not products:
        raise SgError(
            "No products were formed.\n"
            "Possible causes: incompatible ends, no enzyme site, "
            "or the overlap is too short for Gibson/In-Fusion (lower it with --limit)."
        )
    return products


def describe_products(products: list[Dseqrecord]) -> str:
    lines = [f"{len(products)} possible products:"]
    for i, p in enumerate(products, 1):
        top = "circular" if p.circular else "linear"
        lines.append(f"  {i}. {fmt_bp(len(p))}  ({top})")
    if len(products) > 1:
        lines.append(
            "\nNote: Small products are usually unwanted by-products "
            "(a fragment closing on itself). The construct you want is usually the largest."
        )
    return "\n".join(lines)


def sort_products(products: list[Dseqrecord]) -> list[Dseqrecord]:
    """Largest to smallest: the desired construct is usually the largest product."""
    return sorted(products, key=len, reverse=True)


# --------------------------------------------------------------------------
# SnapGene cloning history
# --------------------------------------------------------------------------


def read_history(path: str | Path):
    """Read the cloning history embedded in a SnapGene .dna file."""
    from pydna.snapgene_history_parser import parse_snapgene_history

    path = Path(path)
    if path.suffix.lower() != ".dna":
        raise SgError("Cloning history is only found in SnapGene .dna files.")
    data = path.read_bytes()
    try:
        return parse_snapgene_history(data, file_name=path.name)
    except Exception as e:
        raise SgError(
            f"Could not read history: {e}\n"
            "This file may not have a recorded cloning history "
            "(e.g. imported from elsewhere or created by hand)."
        ) from e


def format_history(rec) -> str:
    """Convert the history tree into readable text."""
    lines = [f"{rec.name}  ({fmt_bp(len(rec))})", ""]
    seen: set[int] = set()

    def walk(node, depth: int = 0) -> None:
        if node is None or id(node) in seen:
            return
        seen.add(id(node))
        pad = "  " * depth
        src = getattr(node, "source", None)
        name = getattr(node, "name", None) or getattr(node, "id", "?")
        size = f" ({fmt_bp(len(node))})" if hasattr(node, "__len__") else ""
        if src is not None:
            kind = type(src).__name__.replace("Source", "")
            lines.append(f"{pad}← {kind}: {name}{size}")
        else:
            lines.append(f"{pad}• {name}{size}")
        for parent in getattr(node, "parent_sequences", None) or []:
            walk(parent, depth + 1)

    walk(rec)
    if len(lines) <= 2:
        lines.append("(history tree appears to be empty)")
    return "\n".join(lines)


# --------------------------------------------------------------------------
# helper: load a list of files
# --------------------------------------------------------------------------


def load_fragments(paths: list[str], linear: bool = False) -> list[Dseqrecord]:
    frags = []
    for p in paths:
        d = load(p, circular=False if linear else None)
        frags.append(d)
    return frags

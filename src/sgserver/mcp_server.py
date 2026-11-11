"""MCP server exposing sgserver to AI agents over stdio.

Run with ``sg-mcp`` (needs the optional dependency: ``pip install "sgserver[mcp]"``).

Safety: every file path an agent passes is resolved and must stay inside the
working root (``SG_MCP_ROOT``, default: the directory the server was started
in). The server only talks over local stdio and opens no network port.
"""

from __future__ import annotations

import functools
import os
from pathlib import Path

from . import api
from .annotate import MIN_IDENTITY
from .core import SgError

try:
    from mcp.server.fastmcp import FastMCP
    from mcp.server.fastmcp.exceptions import ToolError
except ImportError as e:  # pragma: no cover
    raise SystemExit(
        "The MCP server needs the optional dependency:\n"
        '  pip install "sgserver[mcp]"   (or: uv tool install "sgserver[mcp]")'
    ) from e

ROOT = Path(os.environ.get("SG_MCP_ROOT", Path.cwd())).resolve()

mcp = FastMCP(
    "sgserver",
    instructions=(
        "Tools for working with DNA sequence files (.dna, .gb, .fasta, .embl): "
        "inspect, digest, clone, design primers, edit, annotate and verify. "
        "Positions are 1-based and inclusive. File paths must be inside the "
        f"working root ({ROOT}); relative paths are resolved against it. "
        "To clone, pass whole plasmids plus the enzyme names — do not write "
        "digested fragments to files first (GenBank loses sticky ends)."
    ),
)


def _path(p: str | None) -> str | None:
    """Resolve a path and make sure it stays inside ROOT."""
    if p is None:
        return None
    full = (ROOT / p).resolve() if not os.path.isabs(p) else Path(p).resolve()
    if full != ROOT and ROOT not in full.parents:
        raise ToolError(f"Path is outside the working root ({ROOT}): {p}")
    return str(full)


def _paths(ps: list[str]) -> list[str]:
    return [_path(p) for p in ps]


def tool(fn):
    """Register fn as an MCP tool; turn SgError into a clean tool error."""

    @functools.wraps(fn)
    def wrapper(*a, **kw):
        try:
            return fn(*a, **kw)
        except SgError as e:
            raise ToolError(str(e)) from e

    return mcp.tool()(wrapper)


@tool
def sg_status() -> dict:
    """Report whether SnapGene and its bundled tools are available."""
    return api.status()


@tool
def sg_info(file: str, limit: int = 40) -> dict:
    """Summarise a sequence file: length, topology, GC% and features."""
    return api.info(_path(file), limit)


@tool
def sg_convert(file: str, out: str) -> dict:
    """Convert a sequence file between .dna/.gb/.fasta/.embl (.dna needs SnapGene)."""
    return api.convert(_path(file), _path(out))


@tool
def sg_map(file: str, out: str, size: str = "1400,1000") -> dict:
    """Render a SnapGene-style plasmid map image (requires SnapGene)."""
    return api.plasmid_map(_path(file), _path(out), size)


@tool
def sg_history(file: str) -> dict:
    """Show the cloning history embedded in a SnapGene .dna file."""
    return api.history(_path(file))


@tool
def sg_sites(file: str, enzymes: list[str] | None = None, max_cuts: int | None = None,
             all_enzymes: bool = False) -> dict:
    """List restriction sites. max_cuts=1 gives unique cutters."""
    return api.sites(_path(file), enzymes, max_cuts, all_enzymes)


@tool
def sg_digest(file: str, enzymes: list[str], out: str | None = None) -> dict:
    """Digest with the given enzymes; returns fragment sizes (largest first) and ends."""
    return api.digest(_path(file), enzymes, _path(out))


@tool
def sg_gel(lanes: list[str], out: str = "gel.png", enzymes: list[str] | None = None) -> dict:
    """Draw a virtual agarose gel. A lane is a file path or 'name=1000,500' band sizes."""
    safe = []
    for spec in lanes:
        name, _, rest = spec.partition("=")
        if _ and not rest.replace(",", "").strip().isdigit():
            safe.append(f"{name}={_path(rest)}")
        elif not _ and not spec.replace(",", "").strip().isdigit():
            safe.append(_path(spec))
        else:
            safe.append(spec)
    return api.gel(safe, enzymes, _path(out))


@tool
def sg_orfs(file: str, min_aa: int = 100, limit: int = 20) -> dict:
    """Find open reading frames on both strands."""
    return api.orfs(_path(file), min_aa, limit)


@tool
def sg_translate(file: str, start: int, end: int, strand: int = 1, table: int = 1) -> dict:
    """Translate a 1-based inclusive region to protein (strand 1 or -1)."""
    return api.translate(_path(file), start, end, strand, table)


@tool
def sg_find(file: str, query: str, limit: int = 30) -> dict:
    """Search a sequence (IUPAC codes allowed) on both strands."""
    return api.find(_path(file), query, limit)


@tool
def sg_clone(files: list[str], method: str, enzymes: list[str] | None = None,
             min_overlap: int = 25, out: str | None = None) -> dict:
    """Assemble fragments. method: restriction, goldengate, gibson, infusion, ligation,
    gateway-bp, gateway-lr. Pass whole plasmids and enzyme names; the largest
    product is usually the wanted construct."""
    return api.clone(_paths(files), method, enzymes, min_overlap, out=_path(out))


@tool
def sg_pcr(file: str, forward: str, reverse: str, min_match: int = 14,
           out: str | None = None) -> dict:
    """Compute the PCR product of a primer pair on a template."""
    return api.pcr_product(_path(file), forward, reverse, min_match, _path(out))


@tool
def sg_design_primers(file: str, start: int, end: int, target_tm: float = 58.0,
                      forward_tail: str = "", reverse_tail: str = "") -> dict:
    """Design a primer pair for a region, optionally with 5' tails (e.g. restriction sites)."""
    return api.design_primers(_path(file), start, end, target_tm, forward_tail, reverse_tail)


@tool
def sg_check_primers(primers: list[str]) -> dict:
    """Quality report (Tm, GC, hairpin, dimers) for one or two primer sequences."""
    return api.check_primers(primers)


@tool
def sg_fold(sequence_or_file: str) -> dict:
    """Secondary structure and free energy (ViennaRNA, requires SnapGene)."""
    seq = sequence_or_file
    if os.path.exists(seq) or "/" in seq or "." in seq:
        seq = _path(seq)
    return api.fold(seq)


@tool
def sg_annotate(file: str, min_identity: float = MIN_IDENTITY, keep_existing: bool = False,
                out: str | None = None) -> dict:
    """Detect common features (promoters, tags, resistance genes, ori...)."""
    return api.annotate(_path(file), min_identity, keep_existing, _path(out))


@tool
def sg_mutate(file: str, action: str, start: int, end: int | None = None,
              new: str | None = None, mutagenesis_primers: bool = False,
              out: str | None = None) -> dict:
    """Edit a sequence. action 'sub' (replace at start with `new`), 'ins' (insert `new`
    after start) or 'del' (delete start..end)."""
    return api.mutate(_path(file), action, start, end, new, mutagenesis_primers, _path(out))


@tool
def sg_set_origin(file: str, position: int, out: str | None = None) -> dict:
    """Move the origin of a circular plasmid."""
    return api.set_origin(_path(file), position, _path(out))


@tool
def sg_flip(file: str, out: str | None = None) -> dict:
    """Reverse-complement the whole sequence."""
    return api.flip(_path(file), _path(out))


@tool
def sg_add_feature(file: str, start: int, end: int, label: str,
                   kind: str = "misc_feature", strand: int = 1,
                   out: str | None = None) -> dict:
    """Add a labelled feature over a region."""
    return api.add_feature(_path(file), start, end, label, kind, strand, _path(out))


@tool
def sg_verify_read(file: str, read: str) -> dict:
    """Align a sequencing read (.ab1, fasta or gb) to a reference; lists substitutions,
    insertions and deletions."""
    return api.verify_read(_path(file), _path(read))


@tool
def sg_methylation(file: str, enzymes: list[str] | None = None,
                   max_cuts: int | None = None) -> dict:
    """Cut sites overlapping Dam/Dcm methylation, with a known-sensitivity status."""
    return api.methylation(_path(file), enzymes, max_cuts)


@tool
def sg_guides(file: str, pam: str = "NGG", length: int = 20, gc_min: float | None = None,
              gc_max: float | None = None, limit: int = 20) -> dict:
    """CRISPR (SpCas9) guide candidates on both strands."""
    return api.guides(_path(file), pam, length, gc_min, gc_max, limit)


@tool
def sg_align(files: list[str], program: str = "clustalo", out: str | None = None) -> dict:
    """Multiple sequence alignment (clustalo or muscle; requires SnapGene's bundled tools)."""
    return api.align(_paths(files), program, _path(out))


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()

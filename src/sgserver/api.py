"""High-level, JSON-serialisable operations.

Every function here takes plain values (paths, strings, numbers) and returns a
plain ``dict``. The command line (``sg --json``) and the MCP server are thin
layers over this module, so the logic lives in exactly one place.

Conventions: positions are 1-based and inclusive; ``strand`` is 1 or -1;
errors are raised as :class:`sgserver.core.SgError`.
"""

from __future__ import annotations

import re
import tempfile
from pathlib import Path

from . import align as align_mod
from . import analyze, annotate as annotate_mod, clone as clone_mod
from . import edit, methyl, pcr, verify
from .core import (
    SgError,
    load,
    load_many,
    save,
    snapgene_available,
    snapgene_convert,
    snapgene_map,
    topology_str,
)
from .core import snapgene_history as _sg_history


def _feature_row(f) -> dict:
    return {
        "type": f.type,
        "label": analyze.label_of(f),
        "start": int(f.location.start) + 1,
        "end": int(f.location.end),
        "strand": f.location.strand or 0,
    }


def _written(path) -> str:
    return str(Path(path).resolve())


# --------------------------------------------------------------------------
# reading / converting
# --------------------------------------------------------------------------


def info(file: str, limit: int = 40) -> dict:
    """Summary of a sequence file: size, topology, GC and features."""
    rec = load(file)
    seq = str(rec.seq)
    feats = [f for f in rec.features if f.type != "source"]
    return {
        "name": rec.name,
        "description": (rec.description or "").strip(),
        "length": len(seq),
        "topology": topology_str(rec),
        "gc_percent": round(analyze.gc_percent(seq), 2),
        "feature_count": len(feats),
        "features": [_feature_row(f) for f in feats[:limit]],
    }


_CONVERT_FORMATS = {
    ".gb": "GenBank - Standard", ".gbk": "GenBank - Standard",
    ".fa": "FASTA", ".fasta": "FASTA", ".embl": "EMBL",
    ".dna": "SnapGene", ".txt": "Plain Text",
}


def convert(file: str, out: str) -> dict:
    """Convert between sequence formats (.dna needs SnapGene)."""
    src, dst = Path(file), Path(out)
    if dst.suffix.lower() == ".dna" or src.suffix.lower() == ".dna":
        fmt = _CONVERT_FORMATS.get(dst.suffix.lower())
        if fmt is None:
            raise SgError(f"SnapGene cannot convert to this extension: {dst.suffix}")
        snapgene_convert(src, dst, fmt)
    else:
        save(load(src), dst)
    return {"file": _written(dst)}


def plasmid_map(file: str, out: str, size: str = "1400,1000",
                transparent: bool = False, watermark: bool = True) -> dict:
    """Render a SnapGene-style plasmid map image (requires SnapGene)."""
    w, h = (int(x) for x in size.split(","))
    path = snapgene_map(Path(file), Path(out), size=(w, h),
                        transparent=transparent, watermark=watermark)
    return {"file": _written(path)}


def history(file: str, out: str | None = None) -> dict:
    """Cloning history recorded in a SnapGene file."""
    if out:
        return {"file": _written(_sg_history(Path(file), Path(out)))}
    first_error = None
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "history.txt"
        try:
            _sg_history(Path(file), tmp)
            return {"text": tmp.read_text().rstrip()}
        except SgError as e:
            first_error = e
    try:
        rec = clone_mod.read_history(file)
        return {"text": clone_mod.format_history(rec)}
    except SgError:
        raise first_error


# --------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------


def sites(file: str, enzymes: list[str] | None = None, max_cuts: int | None = None,
          all_enzymes: bool = False) -> dict:
    """Restriction enzyme cut positions (max_cuts=1 → unique cutters)."""
    rec = load(file)
    rows = analyze.analyze_sites(rec, enzymes=enzymes, max_cuts=max_cuts,
                                 common_only=not all_enzymes)
    return {
        "length": len(rec.seq),
        "enzymes": [{"name": n, "cuts": len(p), "positions": p} for n, p in rows],
    }


def digest(file: str, enzymes: list[str], out: str | None = None) -> dict:
    """Digest with the given enzymes; fragments are sorted largest first."""
    rec = load(file)
    frags = sorted(analyze.digest(rec, enzymes), key=analyze.fragment_size, reverse=True)
    rows, files = [], []
    for i, f in enumerate(frags, 1):
        row = {"index": i, "size_bp": analyze.fragment_size(f)}
        try:
            k5, o5 = f.seq.five_prime_end()
            k3, o3 = f.seq.three_prime_end()
            row["end_5"] = analyze._end_label(k5, o5)
            row["end_3"] = analyze._end_label(k3, o3)
        except Exception:
            pass
        rows.append(row)
        if out:
            stem = Path(out).with_suffix("")
            files.append(_written(save(f, f"{stem}_{i}_{row['size_bp']}bp.gb")))
    return {"enzymes": enzymes, "fragments": rows, "files": files}


def gel(lanes: list[str], enzymes: list[str] | None = None, out: str = "gel.png") -> dict:
    """Virtual agarose gel. A lane is a file, ``name=file`` or ``name=1000,500``."""
    parsed: dict[str, list[int]] = {}
    for spec in lanes:
        name, rest = spec.split("=", 1) if "=" in spec else (Path(spec).stem, spec)
        if Path(rest).exists():
            rec = load(rest)
            if enzymes:
                parsed[name] = sorted(
                    (analyze.fragment_size(f) for f in analyze.digest(rec, enzymes)),
                    reverse=True)
            else:
                parsed[name] = [len(rec)]
        else:
            parsed[name] = [int(x) for x in rest.split(",")]
    path = analyze.gel_image(parsed, Path(out))
    return {"file": _written(path), "lanes": parsed}


def orfs(file: str, min_aa: int = 100, limit: int = 20) -> dict:
    """Open reading frames on both strands."""
    rec = load(file)
    found = analyze.find_orfs(rec, min_aa=min_aa)
    return {
        "count": len(found),
        "orfs": [
            {"start": s, "end": e, "strand": st, "amino_acids": aa, "length_bp": e - s + 1}
            for s, e, st, aa in found[:limit]
        ],
    }


def translate(file: str, start: int, end: int, strand: int = 1, table: int = 1) -> dict:
    """Translate a 1-based inclusive region to protein."""
    rec = load(file)
    return {
        "name": f"{rec.name}_{start}_{end}",
        "protein": analyze.translate_region(rec, start, end, strand, table),
    }


def find(file: str, query: str, limit: int = 30) -> dict:
    """Search a sequence (IUPAC codes allowed) on both strands."""
    rec = load(file)
    hits = analyze.find_sequence(rec, query)
    return {
        "count": len(hits),
        "hits": [{"start": s, "end": e, "strand": st} for s, e, st in hits[:limit]],
    }


# --------------------------------------------------------------------------
# cloning / PCR / primers
# --------------------------------------------------------------------------


def clone(files: list[str], method: str, enzymes: list[str] | None = None,
          min_overlap: int = 25, linear: bool = False, allow_linear: bool = False,
          out: str | None = None) -> dict:
    """Assemble fragments (restriction, goldengate, gibson, infusion, ligation, gateway)."""
    frags = clone_mod.load_fragments(files, linear=linear)
    products = clone_mod.sort_products(
        clone_mod.assemble(frags, method, enzymes=enzymes, limit=min_overlap,
                           circular_only=not allow_linear)
    )
    written = []
    if out:
        stem = Path(out).stem.split("_")[0] or "construct"
        for p in products:
            p.name = stem[:16]
            p.id = p.name
        if len(products) == 1:
            written.append(_written(save(products[0], out)))
        else:
            base, ext = Path(out).with_suffix(""), Path(out).suffix or ".gb"
            for i, p in enumerate(products, 1):
                written.append(_written(save(p, f"{base}_{i}{ext}")))
    return {
        "method": method,
        "products": [
            {"index": i, "length": len(p),
             "topology": "circular" if p.circular else "linear"}
            for i, p in enumerate(products, 1)
        ],
        "files": written,
    }


def pcr_product(file: str, fwd: str, rev: str, min_match: int = 14,
                out: str | None = None) -> dict:
    """In-silico PCR with a primer pair."""
    rec = load(file)
    amp = pcr.amplify(rec, fwd, rev, limit=min_match)
    f, r = str(amp.forward_primer.seq), str(amp.reverse_primer.seq)
    result = {
        "length": len(amp),
        "sequence": str(amp.seq).upper(),
        "forward_primer": {"sequence": f, "tm": round(pcr.tm(f), 1)},
        "reverse_primer": {"sequence": r, "tm": round(pcr.tm(r), 1)},
        "primer_dimer": pcr.heterodimer(f, r),
    }
    if out:
        result["file"] = _written(save(amp, out))
    return result


def design_primers(file: str, start: int, end: int, target_tm: float = 58.0,
                   fwd_tail: str = "", rev_tail: str = "") -> dict:
    """Design a primer pair for a region (optionally with 5' tails)."""
    rec = load(file)
    d = pcr.design_primers(rec, start, end, target_tm=target_tm,
                           fwd_tail=fwd_tail, rev_tail=rev_tail)
    return {
        "region": {"start": d["region"][0], "end": d["region"][1]},
        "product_bp": d["product_bp"],
        "forward": pcr.primer_report(d["fwd"]),
        "reverse": pcr.primer_report(d["rev"]),
        "primer_dimer": pcr.heterodimer(d["fwd"], d["rev"]),
    }


def check_primers(primers: list[str]) -> dict:
    """Quality report (Tm, GC, hairpin, dimers) for one or two primers."""
    result = {"primers": [pcr.primer_report(p) for p in primers]}
    if len(primers) == 2:
        result["primer_dimer"] = pcr.heterodimer(*primers)
    return result


def fold(seq_or_file: str) -> dict:
    """RNA/DNA secondary structure and free energy (requires SnapGene's ViennaRNA)."""
    seq = seq_or_file
    if Path(seq).exists():
        seq = str(load(seq).seq)
    structure, dg = pcr.fold(seq)
    return {"sequence": seq.upper(), "structure": structure, "dg": dg}


# --------------------------------------------------------------------------
# editing
# --------------------------------------------------------------------------


def annotate(file: str, min_identity: float = annotate_mod.MIN_IDENTITY,
             keep: bool = False, out: str | None = None) -> dict:
    """Detect common features (promoters, ORFs, tags, resistance genes...)."""
    rec = load(file)
    n = annotate_mod.annotate(rec, min_identity=min_identity, replace=not keep)
    rows = []
    for f in rec.features:
        if f.type == "source":
            continue
        row = _feature_row(f)
        ident = 100.0
        for note in f.qualifiers.get("note", []):
            m = re.match(r"identity ([\d.]+)%", note)
            if m:
                ident = float(m.group(1))
        row["identity_pct"] = ident
        rows.append(row)
    result = {"detected": n, "features": rows}
    if out:
        result["file"] = _written(save(rec, out))
    return result


def mutate(file: str, action: str, start: int, end: int | None = None,
           new: str | None = None, primers: bool = False,
           out: str | None = None) -> dict:
    """Point mutation / insertion / deletion. action: 'sub' | 'ins' | 'del'."""
    rec = load(file)
    if action == "sub":
        if not new:
            raise SgError("'sub' requires `new`.")
        result = edit.substitute(rec, start, new)
    elif action == "ins":
        if not new:
            raise SgError("'ins' requires `new`.")
        result = edit.insert(rec, start, new)
    elif action == "del":
        if end is None:
            raise SgError("'del' requires `end`.")
        result = edit.delete(rec, start, end)
    else:
        raise SgError(f"Unknown action: {action} (use sub, ins or del)")
    data = {"action": action, "new_length": len(result)}
    if primers and action == "sub":
        data["mutagenesis_primers"] = edit.mutagenesis_primers(rec, start, new)
    if out:
        data["file"] = _written(save(result, out))
    return data


def set_origin(file: str, position: int, out: str | None = None) -> dict:
    """Move the origin of a circular plasmid."""
    result = edit.set_origin(load(file), position)
    data = {"origin": position, "length": len(result)}
    if out:
        data["file"] = _written(save(result, out))
    return data


def flip(file: str, out: str | None = None) -> dict:
    """Reverse-complement the whole sequence."""
    result = edit.flip(load(file))
    data = {"length": len(result)}
    if out:
        data["file"] = _written(save(result, out))
    return data


def add_feature(file: str, start: int, end: int, label: str,
                kind: str = "misc_feature", strand: int = 1,
                out: str | None = None) -> dict:
    """Add a labelled feature over a region."""
    rec = edit.add_feature(load(file), start, end, label, kind, strand)
    data = {"label": label, "start": start, "end": end}
    if out:
        data["file"] = _written(save(rec, out))
    return data


# --------------------------------------------------------------------------
# validation
# --------------------------------------------------------------------------


def verify_read(file: str, read: str) -> dict:
    """Align a sequencing read (.ab1, fasta, gb) to a reference and list differences."""
    return verify.align_read(load(file), verify.read_trace(read))


def methylation(file: str, enzymes: list[str] | None = None,
                max_cuts: int | None = None, all_enzymes: bool = False) -> dict:
    """Cut sites that overlap Dam/Dcm methylation."""
    from Bio.Restriction import AllEnzymes

    rec = load(file)
    rows = analyze.analyze_sites(rec, enzymes=enzymes, max_cuts=max_cuts,
                                 common_only=not all_enzymes)
    by_enzyme = {n: p for n, p in rows}
    lengths = {}
    for n in by_enzyme:
        enz = next((e for e in AllEnzymes if str(e) == n), None)
        lengths[n] = len(enz.site) if enz else 6
    blocked = methyl.blocked_sites(rec, by_enzyme, lengths)
    out = []
    for name in sorted(blocked):
        for cut, kinds in sorted(blocked[name].items()):
            out.append({"enzyme": name, "position": cut,
                        "methylases": sorted(kinds),
                        "status": methyl.classify(name, kinds)})
    return {"overlapping_sites": out}


def guides(file: str, pam: str = "NGG", length: int = 20, gc_min: float | None = None,
           gc_max: float | None = None, limit: int = 20) -> dict:
    """CRISPR (SpCas9) guide candidates on both strands."""
    rec = load(file)
    found = edit.find_guides(rec, pam=pam, length=length)
    if gc_min or gc_max:
        found = [g for g in found if (gc_min or 0) <= g["gc"] <= (gc_max or 100)]
    return {"pam": pam, "count": len(found), "guides": found[:limit]}


def align(files: list[str], program: str = "clustalo", out: str | None = None) -> dict:
    """Multiple sequence alignment (requires SnapGene's bundled aligners)."""
    records = []
    for f in files:
        records.extend(load_many(f))
    text = align_mod.align(records, program=program, out=out)
    data = {"sequences": len(records), "program": program}
    if out:
        data["file"] = _written(out)
    else:
        data["alignment"] = text
    return data


def status() -> dict:
    """What is available on this machine."""
    from .core import SNAPGENE_DIR, TOOLS

    return {
        "snapgene_available": snapgene_available(),
        "snapgene_dir": str(SNAPGENE_DIR),
        "bundled_tools": TOOLS.exists(),
    }


def doctor() -> dict:
    """Check the installation and report which commands are usable.

    Each check is {name, status: ok|warn|fail, detail, fix}. ``ok`` overall
    means the SnapGene-free core works; ``capabilities`` says which command
    groups are available on this machine.
    """
    import importlib.metadata as md
    import os
    import platform
    import sys

    from .core import SNAPGENE_BIN, SNAPGENE_DIR, TOOLS, gui_running

    checks: list[dict] = []

    def add(name, status, detail, fix=""):
        checks.append({"name": name, "status": status, "detail": detail, "fix": fix})

    add("platform", "ok", f"{platform.system()} {platform.machine()}, Python {platform.python_version()}")
    if sys.version_info < (3, 12):
        add("python", "fail", f"Python {platform.python_version()} (needs >= 3.12)",
            "Install Python 3.12+ (uv python install 3.12)")

    for pkg in ("biopython", "pydna", "primer3-py", "matplotlib"):
        try:
            add(f"dep:{pkg}", "ok", md.version(pkg))
        except md.PackageNotFoundError:
            add(f"dep:{pkg}", "fail", "not installed", f"pip install {pkg}")
    try:
        add("dep:mcp", "ok", md.version("mcp"))
    except md.PackageNotFoundError:
        add("dep:mcp", "warn", "not installed (only needed for sg-mcp)",
            'pip install "sgserver[mcp]"')

    sg_dir_ok = SNAPGENE_DIR.exists()
    add("snapgene:dir", "ok" if sg_dir_ok else "warn",
        str(SNAPGENE_DIR) if sg_dir_ok else f"not found: {SNAPGENE_DIR}",
        "" if sg_dir_ok else "Install SnapGene or set SNAPGENE_DIR (optional)")
    cli_ok = SNAPGENE_BIN.exists()
    if sg_dir_ok:
        add("snapgene:cli", "ok" if cli_ok else "warn",
            str(SNAPGENE_BIN) if cli_ok else f"launcher not found: {SNAPGENE_BIN}",
            "" if cli_ok else "Check SNAPGENE_DIR points at the SnapGene install folder")
    lic = Path.home() / ".snapgene" / "snapgene-license.lic"
    lic_ok = lic.exists()
    if sg_dir_ok:
        add("snapgene:license", "ok" if lic_ok else "warn",
            str(lic) if lic_ok else f"missing: {lic}",
            "" if lic_ok else "Activate SnapGene once, or copy your licence file there")
    tools_ok = TOOLS.exists()
    if sg_dir_ok:
        add("snapgene:tools", "ok" if tools_ok else "warn",
            ", ".join(sorted(p.name for p in TOOLS.iterdir())) if tools_ok
            else f"missing: {TOOLS}")
    ftrs = SNAPGENE_DIR / "resources" / "standardCommonFeatures.ftrs"
    ftrs_ok = ftrs.exists()
    if sg_dir_ok:
        add("snapgene:feature-db", "ok" if ftrs_ok else "warn",
            str(ftrs) if ftrs_ok else f"missing: {ftrs}")
    if cli_ok and gui_running():
        add("snapgene:window", "warn", "the SnapGene window is open",
            "Close it: while open, convert/map/history may silently do nothing")

    add("mcp:root", "ok", os.environ.get("SG_MCP_ROOT", f"{Path.cwd()} (default: current directory)"))

    core_ok = not any(c["status"] == "fail" for c in checks)
    snap = cli_ok and lic_ok
    capabilities = {
        "core (info, sites, digest, clone, pcr, primers, mutate, verify, guides, methyl, gel)": core_ok,
        "write/convert .dna, map, history": snap,
        "annotate (feature database)": ftrs_ok,
        "fold, align (bundled tools)": tools_ok,
        "mcp server": any(c["name"] == "dep:mcp" and c["status"] == "ok" for c in checks),
    }
    return {
        "ok": core_ok,
        "checks": checks,
        "capabilities": capabilities,
        "next_steps": [c["fix"] for c in checks if c["fix"]],
    }

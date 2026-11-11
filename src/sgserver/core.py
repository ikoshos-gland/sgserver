"""Shared infrastructure: file reading/writing, SnapGene installation, bundled tools."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pydna.dseqrecord import Dseqrecord

def _default_snapgene_dir() -> Path:
    if sys.platform == "darwin":
        return Path("/Applications/SnapGene.app/Contents")
    if sys.platform == "win32":
        prog = os.environ.get("ProgramFiles", r"C:\Program Files")
        return Path(prog) / "GSL Biotech" / "SnapGene"
    return Path("/opt/gslbiotech/snapgene")


def _snapgene_bin(base: Path) -> Path:
    for rel in ("snapgene.sh", "MacOS/SnapGene", "snapgene.exe", "snapgene"):
        cand = base / rel
        if cand.exists():
            return cand
    if sys.platform == "darwin":
        return base / "MacOS" / "SnapGene"
    if sys.platform == "win32":
        return base / "snapgene.exe"
    return base / "snapgene.sh"


SNAPGENE_DIR = Path(os.environ.get("SNAPGENE_DIR") or _default_snapgene_dir())
SNAPGENE_BIN = _snapgene_bin(SNAPGENE_DIR)
TOOLS = (SNAPGENE_DIR / "Resources" / "Tools"
         if (SNAPGENE_DIR / "Resources" / "Tools").exists()
         else SNAPGENE_DIR / "resources" / "Tools")

EXT_FORMAT = {
    ".dna": "snapgene",
    ".gb": "genbank",
    ".gbk": "genbank",
    ".genbank": "genbank",
    ".fa": "fasta",
    ".fasta": "fasta",
    ".fsa": "fasta",
    ".seq": "fasta",
    ".embl": "embl",
}


class SgError(Exception):
    """Expected error that is shown to the user."""


# --------------------------------------------------------------------------
# file reading / writing
# --------------------------------------------------------------------------


def detect_format(path: Path) -> str:
    fmt = EXT_FORMAT.get(path.suffix.lower())
    if fmt is None:
        raise SgError(
            f"Unrecognized file extension: {path.suffix or '(none)'}\n"
            f"Supported: {', '.join(sorted(EXT_FORMAT))}"
        )
    return fmt


def load(path: str | Path, circular: bool | None = None) -> Dseqrecord:
    """Read a single-record file as a Dseqrecord."""
    path = Path(path)
    if not path.exists():
        raise SgError(f"File not found: {path}")
    fmt = detect_format(path)
    records = list(SeqIO.parse(str(path), fmt))
    if not records:
        raise SgError(f"No sequences in file: {path}")
    if len(records) > 1:
        print(
            f"[note] {path.name} contains {len(records)} sequences, using the first one.",
            file=sys.stderr,
        )
    rec = records[0]
    if circular is None:
        circular = rec.annotations.get("topology", "linear") == "circular"
    d = Dseqrecord(rec, circular=circular)
    placeholder = ("<unknown name>", "<unknown id>", "<unknown description>", "", None)
    # In SnapGene the file name is the plasmid name; the "name" field inside
    # .dna is often just the first word of the description, so it is unreliable.
    if fmt == "snapgene" or rec.name in placeholder:
        d.name = path.stem
    if d.id in placeholder or fmt == "snapgene":
        d.id = path.stem
    if (rec.description or "") in placeholder:
        d.description = ""
    return d


def load_many(path: str | Path) -> list[SeqRecord]:
    """Multi-record file (e.g. alignment input)."""
    path = Path(path)
    if not path.exists():
        raise SgError(f"File not found: {path}")
    records = list(SeqIO.parse(str(path), detect_format(path)))
    if not records:
        raise SgError(f"No sequences in file: {path}")
    placeholder = ("<unknown name>", "<unknown id>", "<unknown description>", "", None)
    for i, r in enumerate(records):
        if r.id in placeholder or r.name in placeholder:
            r.id = path.stem if len(records) == 1 else f"{path.stem}_{i + 1}"
            r.name = r.id
        r.description = ""
    return records


def save(rec: Dseqrecord | SeqRecord, path: str | Path) -> Path:
    """Write as GenBank or FASTA. If .dna is requested, convert with the SnapGene CLI."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = path.suffix.lower()

    # Only SnapGene's own CLI can write its native format: write .gb first, then convert.
    if suffix == ".dna":
        tmp = path.with_suffix(".tmp.gb")
        _write_bio(rec, tmp, "genbank")
        try:
            snapgene_convert(tmp, path, "SnapGene")
        finally:
            tmp.unlink(missing_ok=True)
        return path

    fmt = detect_format(path)
    _write_bio(rec, path, fmt)
    return path


def _write_bio(rec, path: Path, fmt: str) -> None:
    if fmt == "genbank":
        # Note: do not slice rec (rec[:]) — on a Dseqrecord this opens the circle.
        # The GenBank LOCUS line breaks with names longer than 16 characters.
        name = (rec.name or "sequence").replace(" ", "_")[:16]
        rec.name = name
        if not rec.id or rec.id == "<unknown id>":
            rec.id = name
        rec.annotations.setdefault("molecule_type", "DNA")
        topology = "circular" if getattr(rec, "circular", False) else "linear"
        rec.annotations["topology"] = topology
    with open(path, "w") as fh:
        SeqIO.write(rec, fh, fmt)


# --------------------------------------------------------------------------
# SnapGene's own command line
# --------------------------------------------------------------------------


def snapgene_available() -> bool:
    return SNAPGENE_BIN.exists()


def gui_running() -> bool:
    """Is the SnapGene window open? While it is open, some CLI commands are silently swallowed."""
    if shutil.which("pgrep") is None:
        return False
    r = subprocess.run(
        ["pgrep", "-f", "/opt/gslbiotech/snapgene/snapgene$"],
        capture_output=True,
        text=True,
    )
    return r.returncode == 0 and bool(r.stdout.strip())


GUI_HINT = (
    "\nNote: the SnapGene window is currently open. While it is open, the command line "
    "sometimes silently does nothing. Close SnapGene and try again."
)


def _run_snapgene(args: list[str], timeout: int = 180) -> str:
    if not snapgene_available():
        raise SgError(f"SnapGene not found: {SNAPGENE_BIN}")
    r = subprocess.run(
        [str(SNAPGENE_BIN), *args],
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    out = (r.stdout + r.stderr).strip()
    if "activate before using" in out:
        raise SgError(
            "The SnapGene license is not active for the CLI.\n"
            "The license file must be at ~/.snapgene/snapgene-license.lic."
        )
    return out


def snapgene_convert(src: Path, dst: Path, fmt: str = "GenBank - Standard") -> Path:
    """Convert a file format with the SnapGene CLI."""
    src, dst = Path(src), Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    out = _run_snapgene(["-c", fmt, "-i", str(src), "-o", str(dst)])
    if not dst.exists():
        raise SgError(f"Conversion failed: {out or '(no output)'}"
                      + (GUI_HINT if gui_running() else ""))
    return dst


def snapgene_map(
    src: Path,
    dst: Path,
    size: tuple[int, int] = (1400, 1000),
    transparent: bool = False,
    watermark: bool = True,
) -> Path:
    """Render a plasmid map image in SnapGene's own style."""
    src, dst = Path(src), Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    args = [
        "--createPreview",
        "-s",
        f"{size[0]},{size[1]}",
        "-i",
        str(src),
        "-o",
        str(dst),
    ]
    if transparent:
        args.append("-t")
    if not watermark:
        args.append("-w")
    out = _run_snapgene(args)
    if not dst.exists():
        raise SgError(f"Could not generate map: {out or '(no output)'}"
                      + (GUI_HINT if gui_running() else ""))
    return dst


def snapgene_history(src: Path, dst: Path, size: tuple[int, int] = (1400, 1000)) -> Path:
    """Export the file's cloning history (.txt or image)."""
    src, dst = Path(src), Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    args = ["--history", "-i", str(src), "-o", str(dst)]
    if dst.suffix.lower() != ".txt":
        args += ["-s", f"{size[0]},{size[1]}"]
    out = _run_snapgene(args)
    if not dst.exists():
        raise SgError(
            f"Could not export history: {out or '(no output)'}\n"
            "This file may not have a cloning history recorded in SnapGene."
            + (GUI_HINT if gui_running() else "")
        )
    return dst


# --------------------------------------------------------------------------
# Helper programs bundled inside SnapGene
# --------------------------------------------------------------------------


def tool_env() -> dict[str, str]:
    """Bundled tools need SnapGene's own libunwind."""
    env = dict(os.environ)
    prev = env.get("LD_LIBRARY_PATH", "")
    env["LD_LIBRARY_PATH"] = f"{SNAPGENE_DIR}:{prev}" if prev else str(SNAPGENE_DIR)
    return env


def run_tool(rel_path: str, args: list[str], timeout: int = 300) -> subprocess.CompletedProcess:
    exe = TOOLS / rel_path
    if not exe.exists():
        raise SgError(f"Bundled tool not found: {exe}")
    return subprocess.run(
        [str(exe), *args],
        capture_output=True,
        text=True,
        env=tool_env(),
        timeout=timeout,
    )


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------


def label_of(feature) -> str:
    q = feature.qualifiers
    for key in ("label", "gene", "product", "note", "standard_name"):
        v = q.get(key)
        if v:
            return v[0] if isinstance(v, list) else str(v)
    return feature.type


def fmt_bp(n: int) -> str:
    return f"{n:,} bp"


def revcomp(s: str) -> str:
    return str(Seq(s).reverse_complement())


def topology_str(rec) -> str:
    return "circular" if getattr(rec, "circular", False) else "linear"

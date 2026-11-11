"""Automatic feature detection — using SnapGene's own common-features database.

The `standardCommonFeatures.ftrs` file that ships with SnapGene contains more
than 1400 known features (promoters, resistance genes, ori, tags, primer
sites...). This module reads that database and scans a sequence, marking
matches as features.

SnapGene does not require exact matches: genes such as AmpR differ by a few
bases from plasmid to plasmid. So a seed-based approximate matching is used
here — a short anchor sequence is searched, then full-length identity is
computed.

The database is read from the user's own SnapGene installation; the parsed
version is cached under ~/.cache for speed and is never copied into the project.
"""

from __future__ import annotations

import gzip
import json
import os
from collections import defaultdict
from pathlib import Path

from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation

from .core import SNAPGENE_DIR, SgError, revcomp

FTRS_PATH = SNAPGENE_DIR / "resources" / "standardCommonFeatures.ftrs"

SEED = 16          # anchor length
MIN_IDENTITY = 0.85  # acceptance threshold (0.80 gave false positives like "bom")
MIN_LEN = 15         # shorter features are skipped (M13 primers are 17 bp)


def _cache_path() -> Path:
    base = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
    return base / "sgserver" / "common_features.json.gz"


def build_database() -> list[dict]:
    """Parse the .ftrs file and produce a [{name, type, seq}] list."""
    if not FTRS_PATH.exists():
        raise SgError(
            f"SnapGene feature database not found: {FTRS_PATH}\n"
            "SnapGene must be installed for automatic detection."
        )
    from pydna.snapgene_history_parser import SgffReader

    obj = SgffReader.from_file(str(FTRS_PATH))
    reference = obj.blocks[0][0]["sequence"].upper()

    out: list[dict] = []
    for f in obj.features:
        seq = "".join(
            reference[seg.start : seg.end] for seg in sorted(f.segments, key=lambda s: s.start)
        )
        if len(seq) < MIN_LEN:
            continue
        out.append({"name": f.name, "type": f.type or "misc_feature", "seq": seq})
    return out


def load_database(refresh: bool = False) -> list[dict]:
    cache = _cache_path()
    if cache.exists() and not refresh:
        try:
            with gzip.open(cache, "rt", encoding="utf-8") as fh:
                return json.load(fh)
        except Exception:
            pass  # corrupt cache: rebuild
    db = build_database()
    cache.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(cache, "wt", encoding="utf-8") as fh:
        json.dump(db, fh)
    return db


def _identity(a: str, b: str) -> float:
    if len(a) != len(b) or not a:
        return 0.0
    same = sum(1 for x, y in zip(a, b) if x == y)
    return same / len(a)


def _seed_index(haystack: str) -> dict[str, list[int]]:
    """Positions of every SEED-length window in the haystack."""
    idx: dict[str, list[int]] = defaultdict(list)
    for i in range(len(haystack) - SEED + 1):
        idx[haystack[i : i + SEED]].append(i)
    return idx


def find_features(
    rec,
    min_identity: float = MIN_IDENTITY,
    db: list[dict] | None = None,
) -> list[dict]:
    """Find recognized features; returns [{name,type,start,end,strand,identity}].

    start/end are 0-based, end exclusive (Biopython convention).
    """
    db = db if db is not None else load_database()
    seq = str(rec.seq).upper()
    n = len(seq)
    if n < SEED:
        return []
    circular = bool(getattr(rec, "circular", False))

    longest = max(len(e["seq"]) for e in db)
    pad = min(longest, n) if circular else 0
    haystack = seq + seq[:pad]
    index = _seed_index(haystack)

    hits: list[dict] = []
    for entry in db:
        s = entry["seq"]
        L = len(s)
        if L < MIN_LEN or L > n:
            continue
        for strand, probe in ((1, s), (-1, revcomp(s))):
            if strand == -1 and probe == s:
                continue  # palindrome: do not count twice
            # try several anchors spread across the sequence
            offsets = sorted({0, max(0, (L - SEED) // 2), max(0, L - SEED)})
            seen_starts: set[int] = set()
            for off in offsets:
                for pos in index.get(probe[off : off + SEED], ()):
                    start = pos - off
                    if start < 0 or start >= n or start + L > len(haystack):
                        continue
                    if start in seen_starts:
                        continue
                    seen_starts.add(start)
                    ident = _identity(probe, haystack[start : start + L])
                    if ident >= min_identity:
                        hits.append(
                            {
                                "name": entry["name"],
                                "type": entry["type"],
                                "start": start,
                                "end": start + L,
                                "strand": strand,
                                "identity": ident,
                            }
                        )
    return _dedupe(hits)


def _dedupe(hits: list[dict]) -> list[dict]:
    """Drop duplicates covering the same region; keep the longest and most similar."""
    # on ties prefer the + strand (SnapGene does the same for half-palindromes)
    hits.sort(key=lambda h: (-(h["end"] - h["start"]), -h["identity"], -h["strand"]))
    kept: list[dict] = []
    for h in hits:
        redundant = False
        for k in kept:
            if h["strand"] != k["strand"]:
                continue
            overlap = min(h["end"], k["end"]) - max(h["start"], k["start"])
            if overlap <= 0:
                continue
            shorter = min(h["end"] - h["start"], k["end"] - k["start"])
            # merge records with the same name or that almost fully overlap
            if h["name"] == k["name"] or overlap >= 0.9 * shorter:
                redundant = True
                break
        if not redundant:
            kept.append(h)
    kept.sort(key=lambda h: (h["start"], -(h["end"] - h["start"])))
    return kept


def annotate(
    rec,
    min_identity: float = MIN_IDENTITY,
    replace: bool = True,
    db: list[dict] | None = None,
) -> int:
    """Scan the sequence, mark recognized features, and return the number added."""
    n = len(rec.seq)
    hits = find_features(rec, min_identity=min_identity, db=db)
    if replace:
        rec.features = [f for f in rec.features if f.type == "source"]
    for h in hits:
        if h["end"] <= n:
            location = SimpleLocation(h["start"], h["end"], strand=h["strand"])
        else:  # feature spanning the origin
            location = CompoundLocation(
                [
                    SimpleLocation(h["start"], n, strand=h["strand"]),
                    SimpleLocation(0, h["end"] - n, strand=h["strand"]),
                ]
            )
        quals = {"label": [h["name"]], "note": ["common feature"]}
        if h["identity"] < 1.0:
            quals["note"].append(f"identity {100 * h['identity']:.1f}%")
        rec.features.append(SeqFeature(location=location, type=h["type"], qualifiers=quals))
    rec.features.sort(key=lambda f: (int(f.location.start), -len(f)))
    return len(hits)

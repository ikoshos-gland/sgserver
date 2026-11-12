"""Shared test infrastructure.

Most tests use SnapGene's own sample files as the "right answer": the sample
project contains both the inputs of an operation and the result SnapGene
produced. Tests therefore do not verify themselves; they compare against an
independent reference. Tests that need those files are skipped when SnapGene
is not installed (set SNAPGENE_DIR if it is installed elsewhere).
"""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio.Seq import Seq

import os

SNAPGENE = Path(os.environ.get("SNAPGENE_DIR", "/opt/gslbiotech/snapgene"))
SAMPLES = SNAPGENE / "resources" / "sampleData" / "Sample project"
DESIGN = SAMPLES / "2. Experimental Design"
ANNOTATE = SAMPLES / "1. Annotate and Visualize" / "2. Annotate"
VALIDATE = SAMPLES / "3. Validating Results"
DNA_FILES = SAMPLES / "4. Sample Files" / "DNA files"

requires_samples = pytest.mark.skipif(
    not SAMPLES.exists(), reason="SnapGene sample files are not installed"
)
requires_snapgene_cli = pytest.mark.skipif(
    not (SNAPGENE / "snapgene.sh").exists(), reason="SnapGene CLI is not installed"
)


def same_circle(a, b) -> bool:
    """Are two circular sequences the same? Origin and strand are ignored."""
    A, B = str(a.seq).upper(), str(b.seq).upper()
    if len(A) != len(B):
        return False
    doubled = A + A
    return B in doubled or str(Seq(B).reverse_complement()) in doubled


def labels(rec) -> dict[str, tuple[int, int]]:
    """{label: (start, end)} — 0-based, end exclusive."""
    out = {}
    for f in rec.features:
        if f.type == "source":
            continue
        v = f.qualifiers.get("label")
        name = (v[0] if isinstance(v, list) else v) if v else f.type
        out[name] = (int(f.location.start), int(f.location.end))
    return out


@pytest.fixture(scope="session")
def feature_db():
    """Common-feature database — loaded once and shared (slow).

    It is read from the SnapGene installation, so it is unavailable when
    SnapGene is not installed (e.g. on CI).
    """
    from sgserver.annotate import FTRS_PATH, load_database

    if not FTRS_PATH.exists():
        pytest.skip("SnapGene feature database is not installed")
    return load_database()


@pytest.fixture
def pIB2():
    if not SAMPLES.exists():
        pytest.skip("SnapGene sample files are not installed")
    from sgserver.core import load

    return load(SAMPLES / "pIB2-SEC13-mEGFP.dna")


@pytest.fixture
def pUC19():
    if not SAMPLES.exists():
        pytest.skip("SnapGene sample files are not installed")
    from sgserver.core import load

    return load(DNA_FILES / "Basic Cloning Vectors" / "pUC19.dna")


requires_bundled_tools = pytest.mark.skipif(
    not (SNAPGENE / "resources" / "Tools").exists(),
    reason="SnapGene's bundled tools (ViennaRNA, clustalo, ...) are not installed",
)

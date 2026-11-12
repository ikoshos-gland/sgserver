"""File reading, writing and format conversion."""

from __future__ import annotations

import pytest

from sgserver.core import SgError, load, load_many, save, snapgene_convert
from tests.conftest import (
    DNA_FILES,
    SAMPLES,
    requires_samples,
    requires_snapgene_cli,
)



def test_snapgene_file_is_read(pIB2):
    assert len(pIB2) == 7235
    assert pIB2.circular
    assert pIB2.name == "pIB2-SEC13-mEGFP"  # from the file name, not the description


def test_features_are_read(pIB2):
    types = {f.type for f in pIB2.features}
    assert {"CDS", "promoter", "rep_origin", "primer_bind"} <= types
    labels = {(f.qualifiers.get("label") or [""])[0] for f in pIB2.features}
    assert "mEGFP" in labels
    assert "AmpR" in labels


def test_circular_flag_is_preserved(tmp_path, pIB2):
    """Regression: slicing a Dseqrecord used to open the circle."""
    out = tmp_path / "rt.gb"
    save(pIB2, out)
    again = load(out)
    assert again.circular, "the saved plasmid came back linear"
    assert str(again.seq).upper() == str(pIB2.seq).upper()


@requires_samples
def test_linear_stays_linear(tmp_path):
    rec = load(DNA_FILES / "Fluorescent Protein Genes" / "mCherry.dna")
    assert not rec.circular
    out = tmp_path / "lin.gb"
    save(rec, out)
    assert not load(out).circular


def test_fasta_round_trip(tmp_path, pIB2):
    out = tmp_path / "x.fasta"
    save(pIB2, out)
    assert str(load(out).seq).upper() == str(pIB2.seq).upper()


def test_long_name_does_not_break_locus_line(tmp_path, pIB2):
    pIB2.name = "a-very-long-plasmid-name-that-does-not-fit-the-locus-line"
    out = tmp_path / "long.gb"
    save(pIB2, out)
    again = load(out)
    assert len(again) == len(pIB2)  # the file is still readable


def test_unknown_extension_raises(tmp_path):
    p = tmp_path / "file.xyz"
    p.write_text(">a\nACGT\n")
    with pytest.raises(SgError, match="Unrecognized file extension"):
        load(p)


def test_missing_file_raises(tmp_path):
    with pytest.raises(SgError, match="not found"):
        load(tmp_path / "missing.gb")


@requires_samples
def test_multi_record_file(tmp_path):
    from Bio.SeqIO import write

    recs = [
        load(DNA_FILES / "Fluorescent Protein Genes" / n)
        for n in ("GFP.dna", "mCherry.dna")
    ]
    p = tmp_path / "multi.fasta"
    write(recs, str(p), "fasta")
    assert len(load_many(p)) == 2


@requires_snapgene_cli
def test_conversion_with_snapgene_cli(tmp_path, pIB2):
    """SnapGene CLI .dna → GenBank; features must be kept."""
    out = tmp_path / "out.gb"
    snapgene_convert(SAMPLES / "pIB2-SEC13-mEGFP.dna", out, "GenBank - Standard")
    again = load(out)
    assert len(again) == 7235
    assert again.circular
    labels = {(f.qualifiers.get("label") or [""])[0] for f in again.features}
    assert "mEGFP" in labels


@requires_snapgene_cli
def test_saving_as_dna(tmp_path, pUC19):
    """save() writes .dna through the SnapGene CLI."""
    out = tmp_path / "test.dna"
    save(pUC19, out)
    assert out.exists()
    again = load(out)
    assert len(again) == len(pUC19)
    assert again.circular

"""Restriction sites, digests, gels, ORFs, translation, search."""

from __future__ import annotations

import pytest

from sgserver import analyze
from sgserver.core import SgError
from tests.conftest import requires_samples



# --------------------------------------------------------------------------
# restriction sites
# --------------------------------------------------------------------------


def test_pUC19_mcs_unique_cutters(pUC19):
    """The pUC19 multiple cloning site sits at the positions known from the literature."""
    rows = dict(analyze.analyze_sites(pUC19, enzymes=["EcoRI", "KpnI", "BamHI",
                                                      "SalI", "HindIII"]))
    assert rows["EcoRI"] == [397]
    assert rows["KpnI"] == [413]
    assert rows["BamHI"] == [418]
    assert rows["SalI"] == [430]
    assert rows["HindIII"] == [448]


def test_unique_cutter_filter(pIB2):
    rows = analyze.analyze_sites(pIB2, max_cuts=1)
    assert all(len(sites) == 1 for _, sites in rows)
    names = {n for n, _ in rows}
    for enz in ("BamHI", "EcoRV", "ApaI", "AflIII", "StuI"):
        assert enz in names


def test_non_cutting_enzyme_is_not_listed(pUC19):
    """An enzyme with no site must not appear in the result."""
    rows = dict(analyze.analyze_sites(pUC19, enzymes=["NotI", "EcoRI"]))
    assert "EcoRI" in rows
    assert "NotI" not in rows  # no NotI site in pUC19


def test_site_across_circular_boundary():
    """On a circular sequence a site spanning the origin must be found."""
    from pydna.dseqrecord import Dseqrecord

    # place GGATCC across the end/start boundary: ...GGA | TCC...
    seq = "TCC" + "A" * 50 + "GGA"
    rec = Dseqrecord(seq, circular=True)
    assert analyze.analyze_sites(rec, enzymes=["BamHI"])
    linear = Dseqrecord(seq, circular=False)
    assert not analyze.analyze_sites(linear, enzymes=["BamHI"])


def test_unknown_enzyme_raises(pUC19):
    with pytest.raises(SgError, match="Unknown enzyme"):
        analyze.analyze_sites(pUC19, enzymes=["NoSuchEnzyme"])


# --------------------------------------------------------------------------
# digest
# --------------------------------------------------------------------------


def test_digest_fragment_sizes(pIB2):
    """Cut positions EcoRV 1139, BamHI 1635 → 496 and 6739 bp."""
    frags = analyze.digest(pIB2, ["BamHI", "EcoRV"])
    sizes = sorted(analyze.fragment_size(f) for f in frags)
    assert sizes == [496, 6739]


def test_fragments_sum_to_plasmid_size(pIB2):
    """Regression: pydna's len() counted sticky ends twice."""
    for enzymes in (["BamHI", "EcoRV"], ["BamHI"], ["EcoRV", "ApaI", "AflIII"]):
        frags = analyze.digest(pIB2, enzymes)
        total = sum(analyze.fragment_size(f) for f in frags)
        assert total == len(pIB2), f"{enzymes}: {total} != {len(pIB2)}"


def test_single_enzyme_opens_the_circle(pIB2):
    """A unique cutter on a circular plasmid → one linear fragment of the same size."""
    frags = analyze.digest(pIB2, ["BamHI"])
    assert len(frags) == 1
    assert analyze.fragment_size(frags[0]) == len(pIB2)


def test_sticky_and_blunt_ends(pIB2):
    """BamHI leaves a sticky end (GATC overhang), EcoRV a blunt end."""
    frags = analyze.digest(pIB2, ["BamHI", "EcoRV"])
    text = analyze.format_digest(frags, ["BamHI", "EcoRV"])
    assert "GATC" in text
    assert "blunt" in text


def test_digest_with_non_cutting_enzyme_raises(pUC19):
    with pytest.raises(SgError, match="does not cut"):
        analyze.digest(pUC19, ["NotI"])


# --------------------------------------------------------------------------
# gel
# --------------------------------------------------------------------------


def test_gel_image_is_produced(tmp_path, pIB2):
    frags = analyze.digest(pIB2, ["BamHI", "EcoRV"])
    out = analyze.gel_image(
        {"digest": [analyze.fragment_size(f) for f in frags]}, tmp_path / "gel.png"
    )
    assert out.exists() and out.stat().st_size > 1000


def test_empty_gel_raises(tmp_path):
    with pytest.raises(SgError, match="No bands"):
        analyze.gel_image({"empty": []}, tmp_path / "g.png")


# --------------------------------------------------------------------------
# ORFs and translation
# --------------------------------------------------------------------------


def test_mEGFP_orf_is_found(pIB2):
    """The mEGFP CDS is 1653..2372, 239 aa — must match the file's annotation."""
    orfs = analyze.find_orfs(pIB2, min_aa=200)
    assert any(start == 1653 and aa == 239 for start, _, _, aa in orfs)


def test_AmpR_is_found_on_the_reverse_strand(pIB2):
    """AmpR 3518..4378 is on the reverse strand; the ORF scan must see it too."""
    orfs = analyze.find_orfs(pIB2, min_aa=200)
    assert any(strand == -1 and aa == 286 for _, _, strand, aa in orfs)


def test_min_aa_filter(pIB2):
    long_ = analyze.find_orfs(pIB2, min_aa=300)
    short = analyze.find_orfs(pIB2, min_aa=100)
    assert len(long_) < len(short)
    assert all(aa >= 300 for *_, aa in long_)


def test_translation_of_GFP_start(pIB2):
    prot = analyze.translate_region(pIB2, 1653, 2372)
    assert prot.startswith("MVSKGEE")   # well-known start of the fluorescent protein family
    assert prot.endswith("*")           # stop codon included
    assert len(prot) == 240             # 239 aa + stop


def test_translation_on_reverse_strand(pIB2):
    """AmpR is on the reverse strand; translation must give beta-lactamase."""
    prot = analyze.translate_region(pIB2, 3518, 4378, strand=-1)
    assert prot.startswith("M")
    assert "*" not in prot[:-1]  # no premature stop in the middle


def test_gc_content(pIB2):
    gc = analyze.gc_percent(str(pIB2.seq))
    assert 47.0 < gc < 48.5


# --------------------------------------------------------------------------
# search
# --------------------------------------------------------------------------


def test_sequence_search(pIB2):
    """BamHI site GGATCC at 1634 — consistent with the digest result."""
    hits = analyze.find_sequence(pIB2, "GGATCC")
    assert (1634, 1639, 1) in hits


def test_palindrome_is_counted_once(pIB2):
    """GGATCC is its own complement; it must not be reported twice."""
    hits = analyze.find_sequence(pIB2, "GGATCC")
    positions = [h[0] for h in hits]
    assert len(positions) == len(set(positions))


def test_search_on_the_reverse_strand(pIB2):
    """A non-palindromic sequence must also be found on the reverse strand."""
    hits = analyze.find_sequence(pIB2, "ATGGTGAGCAAG")  # start of mEGFP, + strand
    assert any(s == 1 for *_, s in hits)
    rc_hits = analyze.find_sequence(pIB2, "CTTGCTCACCAT")  # its reverse complement
    assert any(s == -1 for *_, s in rc_hits)


def test_iupac_ambiguity_code(pIB2):
    """N matches any base: GGATCC results ⊂ GGATNN results."""
    exact = {h[0] for h in analyze.find_sequence(pIB2, "GGATCC")}
    ambiguous = {h[0] for h in analyze.find_sequence(pIB2, "GGATNN")}
    assert exact <= ambiguous


def test_invalid_letter_raises(pIB2):
    with pytest.raises(SgError, match="Invalid letter"):
        analyze.find_sequence(pIB2, "GGATZZ")


def test_search_across_circular_boundary():
    from pydna.dseqrecord import Dseqrecord

    rec = Dseqrecord("CCC" + "A" * 40 + "GGG", circular=True)
    assert analyze.find_sequence(rec, "GGGCCC")  # sequence spanning the origin

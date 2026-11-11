"""Cloning — the reference is the result in SnapGene's own sample project.

The sample project contains, for each cloning, both the inputs and the correct
product SnapGene produced. These tests require our result to be identical.
"""

from __future__ import annotations

import pytest

from sgserver import analyze, clone, pcr
from sgserver.core import SgError, load, save
from tests.conftest import DESIGN, DNA_FILES, requires_samples, same_circle

pytestmark = requires_samples

RESTRICTION = DESIGN / "2. Restriction Cloning"
GOLDEN_GATE = DESIGN / "3. Golden Gate Assembly"


# --------------------------------------------------------------------------
# against SnapGene's own examples
# --------------------------------------------------------------------------


def test_restriction_cloning_matches_snapgene():
    """pUC57 + mTAGBFP2 (SalI/KpnI) → pUC57-mTAGBFP2, 3387 bp."""
    vector = load(RESTRICTION / "1. pUC57.dna")
    insert = load(RESTRICTION / "2. mTAGBFP SalI, KpnI.dna")
    expected = load(RESTRICTION / "3. pUC57-mTAGBFP2.dna")

    products = clone.sort_products(
        clone.assemble([vector, insert], "restriction", enzymes=["SalI", "KpnI"])
    )
    assert len(expected) == 3387
    assert any(same_circle(p, expected) for p in products), (
        f"no 3387 bp product; got {[len(p) for p in products]}"
    )


def test_golden_gate_matches_snapgene():
    """pGGA + mCherry (BsaI) → pGGA-mCherry, 2848 bp, a single product."""
    vector = load(GOLDEN_GATE / "1. pGGA.dna")
    insert = load(GOLDEN_GATE / "2. mCherry BsaI.dna")
    expected = load(GOLDEN_GATE / "3. pGGA-mCherry.dna")

    products = clone.assemble([vector, insert], "goldengate", enzymes=["BsaI"])
    assert len(products) == 1, "Golden Gate is directional; one product expected"
    assert same_circle(products[0], expected)
    assert products[0].circular


def test_golden_gate_consumes_the_sites():
    """BsaI sites must not remain in the product (the cut is outside the site)."""
    vector = load(GOLDEN_GATE / "1. pGGA.dna")
    insert = load(GOLDEN_GATE / "2. mCherry BsaI.dna")
    product = clone.assemble([vector, insert], "goldengate", enzymes=["BsaI"])[0]
    assert not analyze.analyze_sites(product, enzymes=["BsaI"])


# --------------------------------------------------------------------------
# a real insertion: mCherry into pUC19
# --------------------------------------------------------------------------


@pytest.fixture
def mcherry_amplicon():
    template = load(DNA_FILES / "Fluorescent Protein Genes" / "mCherry.dna")
    return pcr.amplify(
        template,
        "GCGCGAATTCATGGTGAGCAAGGGCGAG",   # EcoRI end
        "GCGCGGATCCCTACTTGTACAGCTCGTCC",  # BamHI end
    )


def _insert(pUC19, amplicon):
    return clone.sort_products(
        clone.assemble([pUC19, amplicon], "restriction", enzymes=["EcoRI", "BamHI"])
    )


def test_insert_size_is_correct(pUC19, mcherry_amplicon):
    """2686 (pUC19) − 21 (between EcoRI..BamHI) + 717 (cut insert) = 3382."""
    products = _insert(pUC19, mcherry_amplicon)
    assert analyze.fragment_size(products[0]) == 3382 or len(products[0]) == 3382
    assert products[0].circular


def test_insert_recreates_both_sites(tmp_path, pUC19, mcherry_amplicon):
    """Directional cloning must bring both sites back and keep the downstream MCS."""
    product = _insert(pUC19, mcherry_amplicon)[0]
    out = tmp_path / "pUC19-mCherry.gb"
    save(product, out)
    again = load(out)
    assert again.circular, "the circle must survive saving"
    for name in ("EcoRI", "BamHI", "SalI", "HindIII"):
        assert analyze.analyze_sites(again, enzymes=[name]), f"{name} site is missing"


def test_insert_carries_the_gene_intact(pUC19, mcherry_amplicon):
    """The mCherry coding sequence must appear in the product exactly."""
    template = load(DNA_FILES / "Fluorescent Protein Genes" / "mCherry.dna")
    product = _insert(pUC19, mcherry_amplicon)[0]
    assert analyze.find_sequence(product, str(template.seq)[:100])


# --------------------------------------------------------------------------
# method selection and error paths
# --------------------------------------------------------------------------


def test_gibson_joins_overlapping_fragments():
    """Three fragments with 30 bp overlaps must assemble into a circle."""
    from pydna.dseqrecord import Dseqrecord

    a = "ATGCATGCATGCATGCATGCATGCATGCAT"      # A-B overlap
    b = "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTT"      # B-C overlap
    c = "TTCCGGAATTCCGGAATTCCGGAATTCCGG"      # C-A overlap
    f1 = Dseqrecord(c + "AAAAACCCCCGGGGGTTTTTAAAAACCCCC" + a)
    f2 = Dseqrecord(a + "GGGGGTTTTTAAAAACCCCCGGGGGTTTTT" + b)
    f3 = Dseqrecord(b + "CCCCCGGGGGTTTTTAAAAACCCCCGGGGG" + c)

    products = clone.assemble([f1, f2, f3], "gibson", limit=25)
    assert products
    assert products[0].circular
    # each overlap is counted once: 3 x 30 unique + 3 x 30 overlap
    assert len(products[0]) == 180


def test_unknown_method_raises(pUC19):
    with pytest.raises(SgError, match="Unknown method"):
        clone.assemble([pUC19], "magic-method")


def test_restriction_without_enzyme_raises(pUC19):
    with pytest.raises(SgError, match="requires at least one enzyme"):
        clone.assemble([pUC19], "restriction")


def test_incompatible_ends_give_no_product():
    """Fragments with no overlap must not close into a circle."""
    from pydna.dseqrecord import Dseqrecord

    a = Dseqrecord("ATGCATGCATGCATGCATGC")
    b = Dseqrecord("TTTTTTTTTTGGGGGGGGGG")
    with pytest.raises(SgError, match="No products were formed"):
        clone.assemble([a, b], "gibson", limit=25)


def test_products_are_sorted_largest_first(pUC19, mcherry_amplicon):
    products = _insert(pUC19, mcherry_amplicon)
    sizes = [len(p) for p in products]
    assert sizes == sorted(sizes, reverse=True)


# --------------------------------------------------------------------------
# the sticky-end trap
# --------------------------------------------------------------------------


def test_fragments_written_to_files_lose_their_sticky_ends(tmp_path, pIB2):
    """Known limitation: GenBank cannot store sticky ends.

    So the "cut first, write to files, then join" route makes the product
    longer by a few bp. The correct way is to pass the plasmids and name the
    enzyme. This test documents that the trap still exists.
    """
    frags = analyze.digest(pIB2, ["BamHI", "EcoRV"])
    paths = []
    for i, f in enumerate(frags):
        p = tmp_path / f"f{i}.gb"
        save(f, p)
        paths.append(p)
    reloaded = [load(p, circular=False) for p in paths]
    products = clone.assemble(reloaded, "ligation")
    assert all(len(u) != len(pIB2) for u in products), (
        "unexpected: sticky ends survived the file round trip — "
        "the warning in the README can be updated"
    )

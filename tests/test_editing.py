"""PCR/primers, editing, CRISPR guides, read verification, feature detection, methylation."""

from __future__ import annotations

import pytest
from Bio.Seq import Seq

from conftest import (ANNOTATE, DESIGN, labels, requires_bundled_tools,
                      requires_samples, same_circle)
from sgserver import annotate, edit, methyl, pcr, verify
from sgserver.core import SgError, load, revcomp



def _ampr(rec):
    return labels(rec)["AmpR"]


# --------------------------------------------------------------------------
# primers / PCR
# --------------------------------------------------------------------------


def test_tm_rises_with_length():
    assert pcr.tm("ATGGTGAGCAAGGGCGAGGAG") > pcr.tm("ATGGTGAGCAAGGG")


def test_primer_report_gc():
    rep = pcr.primer_report("GGGGCCCC")
    assert rep["gc"] == 100.0 and rep["length"] == 8


def test_empty_primer_raises():
    with pytest.raises(SgError):
        pcr.primer_report("  ")


def test_designed_primers_amplify_the_region(pUC19):
    d = pcr.design_primers(pUC19, 146, 468)
    amp = pcr.amplify(pUC19, d["fwd"], d["rev"])
    assert len(amp) == 468 - 146 + 1 == d["product_bp"]
    assert str(amp.seq).upper() == str(pUC19.seq)[145:468].upper()


@requires_samples
def test_gfp_pcr_matches_snapgene_amplicon_length():
    """Full-length PCR on the GFP template (SnapGene's PCR sample project)."""
    template = load(DESIGN / "1. PCR" / "1. GFP.dna")
    amplified = load(DESIGN / "1. PCR" / "2. GFP Amplified.dna")
    d = pcr.design_primers(template, 1, len(template))
    amp = pcr.amplify(template, d["fwd"], d["rev"])
    assert len(amp) == len(template)
    assert len(amplified) >= len(template)


def test_primer_tail_is_added_to_the_product(pUC19):
    d = pcr.design_primers(pUC19, 146, 468, fwd_tail="GGATCC", rev_tail="GAATTC")
    amp = pcr.amplify(pUC19, d["fwd"], d["rev"])
    s = str(amp.seq).upper()
    assert s.startswith("GGATCC") and s.endswith(revcomp("GAATTC"))
    assert len(amp) == d["product_bp"] + 12


def test_invalid_primer_range_raises(pUC19):
    with pytest.raises(SgError):
        pcr.design_primers(pUC19, 500, 100)


def test_non_annealing_primer_raises(pUC19):
    with pytest.raises(SgError):
        pcr.amplify(pUC19, "A" * 25, "C" * 25)


@requires_bundled_tools
def test_folding_gives_a_negative_energy():
    struct, dg = pcr.fold("GGGAAACCCTTTGGGAAACCC")
    assert len(struct) == 21 and set(struct) <= set("().")
    assert dg < 0


# --------------------------------------------------------------------------
# editing
# --------------------------------------------------------------------------


def test_point_mutation_keeps_features(pUC19):
    s, _ = _ampr(pUC19)
    out = edit.substitute(pUC19, s + 11, "A")
    assert len(out) == len(pUC19)
    assert labels(out) == labels(pUC19)
    assert str(out.seq)[s + 10].upper() == "A"


def test_substitute_changes_only_the_target(pUC19):
    out = edit.substitute(pUC19, 100, "NNNN")
    a, b = str(pUC19.seq).upper(), str(out.seq).upper()
    assert a[:99] == b[:99] and a[103:] == b[103:] and b[99:103] == "NNNN"


def test_deletion_shifts_downstream_features(pUC19):
    s, e = _ampr(pUC19)
    out = edit.delete(pUC19, 1, 10)
    assert len(out) == len(pUC19) - 10
    assert _ampr(out) == (s - 10, e - 10)


def test_insertion_shifts_downstream_features(pUC19):
    s, e = _ampr(pUC19)
    out = edit.insert(pUC19, 0, "AAAA")
    assert len(out) == len(pUC19) + 4
    assert _ampr(out) == (s + 4, e + 4)
    assert str(out.seq).upper().startswith("AAAA")


def test_insertion_inside_a_feature_drops_it(pUC19):
    s, e = _ampr(pUC19)
    out = edit.insert(pUC19, (s + e) // 2, "ACGT")
    assert "AmpR" not in labels(out)


def test_set_origin_rotates_the_sequence(pUC19):
    out = edit.set_origin(pUC19, 100)
    assert str(out.seq)[:10].upper() == str(pUC19.seq)[99:109].upper()
    assert same_circle(out, pUC19)


def test_origin_on_linear_raises():
    from pydna.dseqrecord import Dseqrecord

    with pytest.raises(SgError):
        edit.set_origin(Dseqrecord("ACGT" * 10, circular=False), 5)


def test_flip_twice_gives_the_same_sequence(pUC19):
    twice = edit.flip(edit.flip(pUC19))
    assert str(twice.seq).upper() == str(pUC19.seq).upper()


def test_flip_is_the_reverse_complement(pUC19):
    assert str(edit.flip(pUC19).seq).upper() == revcomp(str(pUC19.seq).upper())


def test_add_feature(pUC19):
    edit.add_feature(pUC19, 10, 20, "test")
    assert labels(pUC19)["test"] == (9, 20)


@pytest.mark.parametrize("fn,args", [
    (edit.substitute, (0, "A")),
    (edit.substitute, (10, "Z")),
    (edit.insert, (10**7, "A")),
    (edit.insert, (1, "XX")),
    (edit.delete, (50, 10)),
    (edit.add_feature, (0, 5, "x")),
])
def test_invalid_input_raises(pUC19, fn, args):
    with pytest.raises(SgError):
        fn(pUC19, *args)


def test_mutagenesis_primers_are_complementary(pUC19):
    p = edit.mutagenesis_primers(pUC19, 500, "GG")
    assert p["rev"] == revcomp(p["fwd"])
    assert len(p["fwd"]) == 15 + 2 + 15
    seq = str(pUC19.seq).upper()
    assert p["fwd"][:15] == seq[484:499]
    assert p["fwd"][15:17] == "GG" and p["old"] == seq[499:501]


# --------------------------------------------------------------------------
# CRISPR guides
# --------------------------------------------------------------------------


def test_guides_match_real_positions(pUC19):
    seq = str(pUC19.seq).upper() * 2  # for guides that span the origin
    guides = edit.find_guides(pUC19)
    assert guides
    for g in guides:
        site = seq[g["start"] - 1 : g["end"]]
        if g["strand"] == 1:
            assert site == g["guide"]
        else:
            assert revcomp(site) == g["guide"]
        assert g["pam"][1:] == "GG" and len(g["guide"]) == 20


def test_guide_limit(pUC19):
    assert len(edit.find_guides(pUC19, limit=3)) == 3


def test_invalid_pam_raises(pUC19):
    with pytest.raises(SgError):
        edit.find_guides(pUC19, pam="XYZ")


# --------------------------------------------------------------------------
# read verification
# --------------------------------------------------------------------------


def test_exact_read_has_no_differences(pUC19):
    read = str(pUC19.seq).upper()[1000:1500]
    r = verify.align_read(pUC19, read)
    assert r["differences"] == [] and r["strand"] == 1
    assert r["start"] == 1001 and r["identity_pct"] == 100.0


def test_reverse_strand_read_is_found(pUC19):
    read = revcomp(str(pUC19.seq).upper()[1000:1500])
    r = verify.align_read(pUC19, read)
    assert r["strand"] == -1 and r["start"] == 1001 and r["differences"] == []


def test_point_substitution_is_found_at_the_right_position(pUC19):
    seq = str(pUC19.seq).upper()
    read = seq[1000:1500]
    alt = "T" if read[200] != "T" else "A"
    mut = read[:200] + alt + read[201:]
    r = verify.align_read(pUC19, mut)
    assert len(r["differences"]) == 1
    d = r["differences"][0]
    assert d["type"] == "substitution" and d["position"] == 1201
    assert d["reference"] == read[200] and d["read"] == alt


def test_deletion_is_found(pUC19):
    read = str(pUC19.seq).upper()[1000:1500]
    r = verify.align_read(pUC19, read[:250] + read[251:])
    assert [d["type"] for d in r["differences"]] == ["deletion"]


def test_insertion_is_found(pUC19):
    read = str(pUC19.seq).upper()[1000:1500]
    r = verify.align_read(pUC19, read[:250] + "G" + read[250:])
    assert [d["type"] for d in r["differences"]] == ["insertion"]


def test_read_spanning_the_circular_boundary(pUC19):
    seq = str(pUC19.seq).upper()
    read = seq[-200:] + seq[:200]
    r = verify.align_read(pUC19, read)
    assert r["differences"] == [] and r["start"] == len(seq) - 199


def test_too_short_read_raises(pUC19):
    with pytest.raises(SgError):
        verify.align_read(pUC19, "ACGT")


def test_scf_gives_a_clear_error(tmp_path):
    p = tmp_path / "x.scf"
    p.write_bytes(b"")
    with pytest.raises(SgError):
        verify.read_trace(p)


def test_missing_file_raises(tmp_path):
    with pytest.raises(SgError):
        verify.read_trace(tmp_path / "yok.ab1")


# --------------------------------------------------------------------------
# automatic feature detection — SnapGene's own annotated pUC19 is the reference
# --------------------------------------------------------------------------


def test_detection_finds_the_same_features_as_snapgene(feature_db):
    plain = load(ANNOTATE / "1. pUC19 Unannotated.dna")
    ref = load(ANNOTATE / "2. pUC19 Annotated.dna")
    found = {h["name"] for h in annotate.find_features(plain, db=feature_db)}
    expected = {n for n in labels(ref)}
    assert expected <= found, f"missing: {expected - found}"


def test_detected_positions_match(feature_db):
    plain = load(ANNOTATE / "1. pUC19 Unannotated.dna")
    ref = labels(load(ANNOTATE / "2. pUC19 Annotated.dna"))
    got = {h["name"]: (h["start"], h["end"])
           for h in annotate.find_features(plain, db=feature_db)}
    for name in ("AmpR", "ori", "AmpR promoter", "lac promoter", "M13 fwd"):
        assert got[name] == ref[name], name


def test_no_false_positives_on_random_sequence(feature_db):
    import random

    random.seed(1)
    rec = Seq("".join(random.choice("ACGT") for _ in range(3000)))
    from pydna.dseqrecord import Dseqrecord

    assert annotate.find_features(Dseqrecord(rec), db=feature_db) == []


def test_annotate_writes_features_to_the_record(feature_db):
    plain = load(ANNOTATE / "1. pUC19 Unannotated.dna")
    n = annotate.annotate(plain, db=feature_db)
    assert n >= 10 and "AmpR" in labels(plain)


# --------------------------------------------------------------------------
# methylation
# --------------------------------------------------------------------------


def test_ClaI_overlaps_dam():
    from pydna.dseqrecord import Dseqrecord

    rec = Dseqrecord("AAAAAAAAATCGATCAAAAAAAAA", circular=False)  # ATCGATC: GATC inside ClaI
    kesim = {"ClaI": [11]}
    blocked = methyl.blocked_sites(rec, kesim, {"ClaI": 6})
    assert "Dam" in blocked["ClaI"][11]


def test_unmethylated_sequence_returns_empty():
    from pydna.dseqrecord import Dseqrecord

    rec = Dseqrecord("AAAAAAAAAAGAATTCAAAAAAAAAA", circular=False)
    assert methyl.blocked_sites(rec, {"EcoRI": [12]}, {"EcoRI": 6}) == {}


def test_BamHI_overlaps_dam_but_is_unaffected():
    assert methyl.classify("BamHI", {"Dam"}) == "unaffected"
    assert methyl.classify("MboI", {"Dam"}) == "blocked"
    assert methyl.classify("DpnI", {"Dam"}) == "required"
    assert methyl.classify("XyzI", {"Dam"}) == "overlap"


def test_dcm_CCWGG_is_found():
    from pydna.dseqrecord import Dseqrecord

    rec = Dseqrecord("AAAAAACCTGGAAAAAA", circular=False)
    got = methyl.blocked_sites(rec, {"X": [9]}, {"X": 5})
    assert "Dcm" in got["X"][9]

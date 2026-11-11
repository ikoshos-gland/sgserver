"""API and CLI tests that need NO SnapGene installation and no proprietary data.

They use a small synthetic plasmid built on the fly, so they run anywhere
(including CI).
"""

from __future__ import annotations

import json
import random

import pytest

from sgserver import api, cli
from sgserver.core import SgError, load, save


SITES = ("GAATTC", "GGATCC", "AAGCTT")  # all palindromic


def _plasmid_seq() -> str:
    """Synthetic plasmid with exactly one EcoRI, BamHI and HindIII site."""
    rng = random.Random(7)

    def filler(n):
        while True:  # regenerate until the filler contains none of the sites
            s = "".join(rng.choice("ACGT") for _ in range(n))
            if not any(site in s for site in SITES):
                return s

    seq = (filler(300) + SITES[0] + filler(120) + SITES[1] + filler(150)
           + SITES[2] + filler(400))
    for site in SITES:
        assert seq.count(site) == 1, site
    return seq


@pytest.fixture
def plasmid(tmp_path):
    from pydna.dseqrecord import Dseqrecord

    rec = Dseqrecord(_plasmid_seq(), circular=True)
    rec.name = rec.id = "synth"
    path = tmp_path / "synth.gb"
    save(rec, path)
    return path


def test_info(plasmid):
    r = api.info(str(plasmid))
    assert r["length"] == len(_plasmid_seq())
    assert r["topology"] == "circular"
    assert 0 < r["gc_percent"] < 100


def test_sites_finds_planted_unique_cutters(plasmid):
    r = api.sites(str(plasmid), enzymes=["EcoRI", "BamHI", "HindIII"])
    names = {e["name"]: e["positions"] for e in r["enzymes"]}
    assert set(names) == {"EcoRI", "BamHI", "HindIII"}
    assert all(e["cuts"] == 1 for e in r["enzymes"])


def test_digest_fragments_sum_to_length(plasmid):
    r = api.digest(str(plasmid), ["EcoRI", "BamHI"])
    assert len(r["fragments"]) == 2
    assert sum(f["size_bp"] for f in r["fragments"]) == len(_plasmid_seq())
    sizes = [f["size_bp"] for f in r["fragments"]]
    assert sizes == sorted(sizes, reverse=True)


def test_find_reports_both_strands(plasmid):
    r = api.find(str(plasmid), "GAATTC")
    assert r["count"] == 1 and r["hits"][0]["strand"] == 1


def test_mutate_substitution(plasmid, tmp_path):
    out = tmp_path / "m.gb"
    r = api.mutate(str(plasmid), "sub", 10, new="NN", out=str(out))
    assert r["new_length"] == len(_plasmid_seq())
    assert str(load(out).seq).upper()[9:11] == "NN"


def test_mutate_requires_new(plasmid):
    with pytest.raises(SgError):
        api.mutate(str(plasmid), "sub", 10)


def test_pcr_product_and_design(plasmid):
    d = api.design_primers(str(plasmid), 50, 400)
    r = api.pcr_product(str(plasmid), d["forward"]["sequence"], d["reverse"]["sequence"])
    assert r["length"] == d["product_bp"] == 351
    assert r["sequence"] == _plasmid_seq()[49:400]


def test_cloning_with_restriction(plasmid, tmp_path):
    """Cut the plasmid and an amplicon with the same enzymes and join them."""
    d = api.design_primers(str(plasmid), 100, 420,
                           fwd_tail="GCGCGAATTC", rev_tail="GCGCGGATCC")
    amp = tmp_path / "amp.gb"
    api.pcr_product(str(plasmid), d["forward"]["sequence"], d["reverse"]["sequence"],
                    out=str(amp))
    r = api.clone([str(plasmid), str(amp)], "restriction", enzymes=["EcoRI", "BamHI"])
    assert r["products"], "no cloning product"
    assert r["products"][0]["topology"] == "circular"


def test_verify_read(plasmid):
    seq = _plasmid_seq()
    read = seq[100:400]
    alt = "A" if read[150] != "A" else "C"
    mutated = read[:150] + alt + read[151:]
    reads = plasmid.parent / "read.fasta"
    reads.write_text(f">r\n{mutated}\n")
    r = api.verify_read(str(plasmid), str(reads))
    assert [d["type"] for d in r["differences"]] == ["substitution"]
    assert r["differences"][0]["position"] == 251


def test_guides_have_pam(plasmid):
    r = api.guides(str(plasmid), limit=5)
    assert r["count"] >= 5 and len(r["guides"]) == 5
    assert all(g["pam"].endswith("GG") for g in r["guides"])


def test_methylation_known_cases(plasmid):
    r = api.methylation(str(plasmid), enzymes=["BamHI"])
    # BamHI overlaps Dam (GGATCC contains GATC) but is unaffected
    assert r["overlapping_sites"][0]["status"] == "unaffected"


# --------------------------------------------------------------------------
# command line
# --------------------------------------------------------------------------


def _run(capsys, *argv):
    code = cli.main(list(argv))
    out = capsys.readouterr()
    return code, out.out, out.err


def test_cli_json_is_valid(capsys, plasmid):
    code, out, _ = _run(capsys, "--json", "sites", str(plasmid), "-e", "EcoRI")
    assert code == 0
    data = json.loads(out)
    assert data["enzymes"][0]["name"] == "EcoRI"


def test_cli_json_error_is_structured(capsys, plasmid):
    code, out, _ = _run(capsys, "--json", "digest", str(plasmid), "-e", "NoSuchEnzyme")
    assert code == 1
    assert "Unknown enzyme" in json.loads(out)["error"]


def test_cli_text_error_goes_to_stderr(capsys, plasmid):
    code, out, err = _run(capsys, "digest", str(plasmid), "-e", "NoSuchEnzyme")
    assert code == 1 and out == "" and "Unknown enzyme" in err


@pytest.mark.parametrize("cmd", [
    ["info"], ["sites"], ["orf", "--min-aa", "20"], ["guides", "--limit", "3"],
    ["methyl"], ["primers", "50", "400"],
])
def test_cli_every_json_command_parses(capsys, plasmid, cmd):
    name, *rest = cmd
    code, out, _ = _run(capsys, "--json", name, str(plasmid), *rest)
    assert code == 0
    json.loads(out)


def test_cli_help_is_english(capsys):
    code, out, _ = _run(capsys)
    assert code == 0 and "Read / convert" in out


def test_status_does_not_require_snapgene(capsys):
    code, out, _ = _run(capsys, "--json", "status")
    assert code == 0
    assert "snapgene_available" in json.loads(out)


def test_doctor_and_conf_alias(capsys):
    code, out, _ = _run(capsys, "--json", "conf")
    assert code == 0
    data = json.loads(out)
    assert data["ok"] is True
    assert "capabilities" in data and "checks" in data

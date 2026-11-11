"""MCP server: tool registration, path sandbox and error handling (no SnapGene needed)."""

from __future__ import annotations

import asyncio

import pytest

pytest.importorskip("mcp")

from tests.test_api import _plasmid_seq  # noqa: E402


@pytest.fixture
def server(tmp_path, monkeypatch):
    monkeypatch.setenv("SG_MCP_ROOT", str(tmp_path))
    import importlib

    from sgserver import mcp_server

    return importlib.reload(mcp_server)


def _call(server, name, args):
    return asyncio.run(server.mcp.call_tool(name, args))


def test_all_tools_are_registered(server):
    tools = asyncio.run(server.mcp.list_tools())
    names = {t.name for t in tools}
    assert len(names) >= 20
    assert {"sg_info", "sg_digest", "sg_clone", "sg_design_primers",
            "sg_verify_read", "sg_status"} <= names


def test_tool_works_inside_root(server, tmp_path):
    from pydna.dseqrecord import Dseqrecord

    from sgserver.core import save

    rec = Dseqrecord(_plasmid_seq(), circular=True)
    rec.name = rec.id = "synth"
    save(rec, tmp_path / "synth.gb")
    out = _call(server, "sg_sites", {"file": "synth.gb", "enzymes": ["EcoRI"]})
    assert "EcoRI" in str(out)


def test_path_outside_root_is_rejected(server):
    with pytest.raises(Exception, match="outside the working root"):
        _call(server, "sg_info", {"file": "/etc/passwd"})
    with pytest.raises(Exception, match="outside the working root"):
        _call(server, "sg_info", {"file": "../../etc/passwd"})


def test_sg_errors_become_tool_errors(server, tmp_path):
    (tmp_path / "x.fasta").write_text(">a\nACGTACGTACGT\n")
    with pytest.raises(Exception, match="Unknown enzyme"):
        _call(server, "sg_digest", {"file": "x.fasta", "enzymes": ["NoSuchEnzyme"]})

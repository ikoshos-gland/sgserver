"""sg — command line toolkit for SnapGene-style sequence work.

Every command accepts the global ``--json`` flag (put it before the command):
results are printed as JSON and errors as ``{"error": "..."}``.
"""

from __future__ import annotations

import argparse
import json
import sys

from . import annotate as annotate_mod
from . import api
from .core import SgError, fmt_bp

HELP = """\
sg — SnapGene toolkit (add --json before the command for machine-readable output)

Read / convert
  info      file summary: length, topology, GC, features
  convert   convert format (.dna <-> .gb/.fasta/.embl)
  map       plasmid map image in SnapGene style
  history   show the cloning history embedded in a SnapGene file

Analysis
  sites     restriction enzyme sites (incl. unique cutters)
  digest    digest and report fragment sizes
  gel       draw a virtual agarose gel
  orf       find open reading frames
  translate translate a region to protein
  find      search a sequence (IUPAC supported)

Cloning / PCR
  clone     assemble fragments (ligation, gibson, goldengate, gateway...)
  pcr       compute the product of a primer pair
  primers   design primers for a region
  check     primer quality: Tm, GC, hairpin, dimer
  fold      secondary structure and free energy (ViennaRNA)

Editing
  annotate  auto-detect common features (SnapGene database)
  mutate    point mutation, insertion, deletion
  origin    move the origin of a circular plasmid
  flip      reverse complement
  feature   add a feature by hand

Validation
  verify    align a sequencing read to a reference, list differences
  methyl    cut sites overlapping Dam/Dcm methylation
  guides    CRISPR guide candidates (SpCas9)

Alignment
  align     multiple alignment (clustalo / muscle)

Setup
  doctor    check the installation and list what works (alias: conf)
  status    short availability summary
"""


def _strand(s: int) -> str:
    return {1: "+", -1: "-"}.get(s, " ")


def _dimer_line(d: dict) -> str:
    if d["tm"] is None:
        return "Cross dimer: none"
    return f"Cross dimer: Tm {d['tm']:.1f} °C, ΔG {d['dg']:.1f} kcal/mol"


def _primer_text(rep: dict, name: str) -> str:
    lines = [
        name,
        f"  Sequence    : {rep['sequence']}",
        f"  Length      : {rep['length']} nt",
        f"  Tm          : {rep['tm']:.1f} °C",
        f"  GC          : {rep['gc']:.1f}%",
        f"  3' GC clamp : {rep['gc_clamp_3p']}/5",
    ]
    for label, key in (("Hairpin     ", "hairpin"), ("Self-dimer  ", "homodimer")):
        if rep[f"{key}_tm"] is None:
            lines.append(f"  {label}: none")
        else:
            warn = "  ⚠ high" if rep[f"{key}_tm"] > 45 else ""
            lines.append(f"  {label}: Tm {rep[f'{key}_tm']:.1f} °C, "
                         f"ΔG {rep[f'{key}_dg']:.1f} kcal/mol{warn}")
    return "\n".join(lines)


def _wrote(r: dict) -> list[str]:
    files = r.get("files") or ([r["file"]] if r.get("file") else [])
    return [f"Written: {f}" for f in files]


# --------------------------------------------------------------------------
# text renderers: result dict -> human-readable text
# --------------------------------------------------------------------------


def _t_info(r, a):
    lines = [
        f"Name      : {r['name']}",
        f"Length    : {fmt_bp(r['length'])}",
        f"Topology  : {r['topology']}",
        f"GC content: {r['gc_percent']:.1f}%",
    ]
    if r["description"] and r["description"] != r["name"]:
        lines.append(f"Description: {r['description']}")
    lines.append(f"\nFeatures ({r['feature_count']}):")
    if not r["features"]:
        lines.append("  (none)")
    for f in r["features"]:
        lines.append(f"  {f['start']:>6}..{f['end']:<6} {_strand(f['strand'])}  "
                     f"{f['type']:<14} {f['label']}")
    if r["feature_count"] > len(r["features"]):
        lines.append(f"  ... and {r['feature_count'] - len(r['features'])} more")
    return "\n".join(lines)


def _t_sites(r, a):
    title = "Unique cutters" if a.max_cuts == 1 else "Restriction sites"
    if not r["enzymes"]:
        return f"{title}: (none)"
    lines = [f"{title} ({len(r['enzymes'])} enzymes):"]
    for e in r["enzymes"]:
        pos = ", ".join(str(s) for s in e["positions"][:12])
        more = f" ... (+{e['cuts'] - 12})" if e["cuts"] > 12 else ""
        lines.append(f"  {e['name']:<12} {e['cuts']:>3} cuts   {pos}{more}")
    return "\n".join(lines)


def _t_digest(r, a):
    lines = [f"{' + '.join(r['enzymes'])} digest → {len(r['fragments'])} fragments:"]
    for f in r["fragments"]:
        ends = f"   ends: {f['end_5']} / {f['end_3']}" if "end_5" in f else ""
        lines.append(f"  {f['index']:>2}. {fmt_bp(f['size_bp']):>12}{ends}")
    lines += [f"Written: {p}" for p in r["files"]]
    return "\n".join(lines)


def _t_orf(r, a):
    if not r["count"]:
        return f"No ORFs above {a.min_aa} aa."
    lines = [f"{r['count']} ORFs (≥{a.min_aa} aa):"]
    for o in r["orfs"]:
        lines.append(f"  {o['start']:>6}..{o['end']:<6} {_strand(o['strand'])}  "
                     f"{o['amino_acids']:>4} aa  ({fmt_bp(o['length_bp'])})")
    return "\n".join(lines)


def _t_find(r, a):
    if not r["count"]:
        return "No matches."
    lines = [f"{r['count']} matches:"]
    lines += [f"  {h['start']:>6}..{h['end']:<6} {_strand(h['strand'])}" for h in r["hits"]]
    return "\n".join(lines)


def _t_clone(r, a):
    ps = r["products"]
    lines = [f"{len(ps)} possible products:"]
    lines += [f"  {p['index']}. {fmt_bp(p['length'])}  ({p['topology']})" for p in ps]
    if len(ps) > 1:
        lines.append("\nNote: small products are usually unwanted by-products "
                     "(a fragment closing on itself). The construct you want "
                     "is usually the largest.")
    return "\n".join(lines + _wrote(r))


def _t_pcr(r, a):
    f, rv = r["forward_primer"], r["reverse_primer"]
    lines = [
        f"PCR product: {fmt_bp(r['length'])}",
        f"  Forward primer: {f['sequence']}  (Tm {f['tm']:.1f} °C)",
        f"  Reverse primer: {rv['sequence']}  (Tm {rv['tm']:.1f} °C)",
    ]
    d = r["primer_dimer"]
    if d["tm"] is None:
        lines.append("  Primer dimer: none")
    else:
        lines.append(f"  Primer dimer: Tm {d['tm']:.1f} °C, ΔG {d['dg']:.1f}"
                     + ("  ⚠" if d["tm"] > 45 else ""))
    return "\n".join(lines + _wrote(r))


def _t_primers(r, a):
    reg = r["region"]
    return "\n".join([
        f"Target: {reg['start']}..{reg['end']}  → product {fmt_bp(r['product_bp'])}\n",
        _primer_text(r["forward"], "Forward primer"), "",
        _primer_text(r["reverse"], "Reverse primer"), "",
        _dimer_line(r["primer_dimer"]),
    ])


def _t_check(r, a):
    parts = []
    for i, p in enumerate(r["primers"], 1):
        parts += [_primer_text(p, f"Primer {i}"), ""]
    if "primer_dimer" in r:
        parts.append(_dimer_line(r["primer_dimer"]))
    return "\n".join(parts).rstrip()


def _t_annotate(r, a):
    lines = [f"{r['detected']} features detected."]
    for f in r["features"]:
        ident = "exact match" if f["identity_pct"] >= 100 else f"identity {f['identity_pct']:.1f}%"
        lines.append(f"  {f['start']:>6}..{f['end']:<6} {_strand(f['strand'])} "
                     f"{f['type']:<14} {f['label']:<22} {ident}")
    return "\n".join(lines + ([""] if r.get("file") else []) + _wrote(r))


def _t_mutate(r, a):
    what = {"sub": f"Substituted {len(a.new or '')} bp at {a.start}",
            "ins": f"Inserted {len(a.new or '')} bp after {a.start}",
            "del": f"Deleted {a.start}..{a.end}"}[r["action"]]
    lines = [what, f"New length: {fmt_bp(r['new_length'])}"]
    p = r.get("mutagenesis_primers")
    if p:
        lines += [f"\nMutagenesis primers ({p['old']} → {p['new']}):",
                  f"  forward: {p['fwd']}", f"  reverse: {p['rev']}"]
    return "\n".join(lines + _wrote(r))


def _t_verify(r, a):
    direction = "forward (+)" if r["strand"] == 1 else "reverse (−)"
    lines = [
        f"Read of {r['read_length']} bp aligned on the {direction} strand",
        f"Reference {r['start']}..{r['end']} "
        f"({r['aligned']} bp overlap, identity {r['identity_pct']:.2f}%)",
    ]
    if not r["differences"]:
        return "\n".join(lines + ["No differences — the read matches the reference exactly."])
    lines.append(f"\n{len(r['differences'])} differences:")
    for d in r["differences"][:25]:
        lines.append(f"  {d['position']:>6}  {d['type']:<12} {d['reference']} → {d['read']}")
    if len(r["differences"]) > 25:
        lines.append(f"  ... and {len(r['differences']) - 25} more")
    return "\n".join(lines)


_METHYL_LABEL = {
    "blocked": "BLOCKED", "required": "methylation REQUIRED",
    "unaffected": "unaffected", "overlap": "overlaps — check the enzyme's sensitivity",
}


def _t_methyl(r, a):
    if not r["overlapping_sites"]:
        return "No cut sites overlap Dam/Dcm sites."
    lines = ["Cut sites overlapping Dam/Dcm methylation:"]
    for s in r["overlapping_sites"]:
        lines.append(f"  {s['enzyme']:<12} {s['position']:>6}   "
                     f"{'+'.join(s['methylases']):<8} {_METHYL_LABEL[s['status']]}")
    lines.append(
        "\nOverlap alone does not mean blocked: GGATCC contains GATC but BamHI is not "
        "affected by Dam.\nCheck the enzyme's sensitivity in REBASE or the supplier "
        "catalogue; if it is sensitive, use a dam-/dcm- strain.")
    return "\n".join(lines)


def _t_guides(r, a):
    lines = [f"{r['count']} guide candidates (PAM {r['pam']}):"]
    for g in r["guides"]:
        lines.append(f"  {g['start']:>6}..{g['end']:<6} {_strand(g['strand'])}  "
                     f"{g['guide']} {g['pam']}  GC {g['gc']:.0f}%")
    return "\n".join(lines)


def _t_align(r, a):
    if r.get("file"):
        return f"Alignment: {r['file']}  ({r['sequences']} sequences)"
    return r["alignment"]


_MARK = {"ok": "✓", "warn": "!", "fail": "✗"}


def _t_doctor(r, a):
    lines = ["sgserver configuration check", ""]
    for c in r["checks"]:
        lines.append(f"  {_MARK[c['status']]} {c['name']:<22} {c['detail']}")
        if c["fix"]:
            lines.append(f"      → {c['fix']}")
    lines += ["", "Available:"]
    for name, ok in r["capabilities"].items():
        lines.append(f"  {'✓' if ok else '✗'} {name}")
    lines.append("")
    lines.append("Result: " + ("core OK" if r["ok"] else "PROBLEMS FOUND — see ✗ above"))
    return "\n".join(lines)


def _t_file(label):
    return lambda r, a: f"{label}: {r['file']}"


# --------------------------------------------------------------------------
# command table: name -> (api call, text renderer)
# --------------------------------------------------------------------------


def _file_cmd(label, fn):
    return (fn, _t_file(label))


COMMANDS = {
    "info": (lambda a: api.info(a.file, a.limit), _t_info),
    "convert": (lambda a: api.convert(a.file, a.out), _t_file("Written")),
    "map": (lambda a: api.plasmid_map(a.file, a.out, a.size, a.transparent,
                                      not a.no_watermark), _t_file("Map")),
    "history": (lambda a: api.history(a.file, a.out),
                lambda r, a: r["text"] if "text" in r else f"History: {r['file']}"),
    "sites": (lambda a: api.sites(a.file, a.enzyme, a.max_cuts, a.all_enzymes), _t_sites),
    "digest": (lambda a: api.digest(a.file, a.enzyme, a.out), _t_digest),
    "gel": (lambda a: api.gel(a.lane, a.enzyme, a.out), _t_file("Gel")),
    "orf": (lambda a: api.orfs(a.file, a.min_aa, a.limit), _t_orf),
    "translate": (lambda a: api.translate(a.file, a.start, a.end, a.strand, a.table),
                  lambda r, a: f">{r['name']}\n{r['protein']}"),
    "find": (lambda a: api.find(a.file, a.query, a.limit), _t_find),
    "clone": (lambda a: api.clone(a.file, a.method, a.enzyme, a.limit, a.linear,
                                  a.allow_linear, a.out), _t_clone),
    "pcr": (lambda a: api.pcr_product(a.file, a.fwd, a.rev, a.limit, a.out), _t_pcr),
    "primers": (lambda a: api.design_primers(a.file, a.start, a.end, a.tm,
                                             a.fwd_tail, a.rev_tail), _t_primers),
    "check": (lambda a: api.check_primers(a.primer), _t_check),
    "fold": (lambda a: api.fold(a.seq),
             lambda r, a: f"{r['sequence']}\n{r['structure']}\nΔG = {r['dg']:.2f} kcal/mol"),
    "annotate": (lambda a: api.annotate(a.file, a.min_identity, a.keep, a.out), _t_annotate),
    "mutate": (lambda a: api.mutate(a.file, a.action, a.start, a.end, a.new,
                                    a.primers, a.out), _t_mutate),
    "origin": (lambda a: api.set_origin(a.file, a.position, a.out),
               lambda r, a: "\n".join([f"Origin moved to position {r['origin']}."] + _wrote(r))),
    "flip": (lambda a: api.flip(a.file, a.out),
             lambda r, a: "\n".join([f"Reverse complement taken ({fmt_bp(r['length'])})."]
                                    + _wrote(r))),
    "feature": (lambda a: api.add_feature(a.file, a.start, a.end, a.label, a.type,
                                          a.strand, a.out),
                lambda r, a: "\n".join([f"Feature added: {r['label']} "
                                        f"({r['start']}..{r['end']})"] + _wrote(r))),
    "verify": (lambda a: api.verify_read(a.file, a.read), _t_verify),
    "methyl": (lambda a: api.methylation(a.file, a.enzyme, a.max_cuts, a.all_enzymes),
               _t_methyl),
    "guides": (lambda a: api.guides(a.file, a.pam, a.length, a.gc_min, a.gc_max, a.limit),
               _t_guides),
    "align": (lambda a: api.align(a.file, a.program, a.out), _t_align),
    "doctor": (lambda a: api.doctor(), _t_doctor),
    "status": (lambda a: api.status(),
               lambda r, a: "\n".join(f"{k}: {v}" for k, v in r.items())),
}


# --------------------------------------------------------------------------
# arguments
# --------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="sg",
        description="SnapGene toolkit",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=HELP,
    )
    p.add_argument("--json", action="store_true",
                   help="print machine-readable JSON (errors as {\"error\": ...})")
    sub = p.add_subparsers(dest="cmd")

    def add(name, help_, aliases=()):
        s = sub.add_parser(name, aliases=list(aliases), help=help_)
        s.set_defaults(cmd=name)
        return s

    s = add("info", "file summary")
    s.add_argument("file")
    s.add_argument("--limit", type=int, default=40, help="how many features to show")

    s = add("convert", "convert format")
    s.add_argument("file")
    s.add_argument("out")

    s = add("map", "plasmid map image (needs SnapGene)")
    s.add_argument("file")
    s.add_argument("out")
    s.add_argument("--size", default="1400,1000")
    s.add_argument("--transparent", action="store_true")
    s.add_argument("--no-watermark", action="store_true")

    s = add("history", "cloning history")
    s.add_argument("file")
    s.add_argument("--out", help="export to a file (.txt or .png)")

    s = add("sites", "restriction sites")
    s.add_argument("file")
    s.add_argument("-e", "--enzyme", action="append", help="specific enzyme(s)")
    s.add_argument("--max-cuts", type=int, help="at most this many cuts (1 = unique cutters)")
    s.add_argument("--all-enzymes", action="store_true", help="all enzymes (slow)")

    s = add("digest", "restriction digest")
    s.add_argument("file")
    s.add_argument("-e", "--enzyme", action="append", required=True)
    s.add_argument("--out", help="write fragments to files")

    s = add("gel", "virtual gel")
    s.add_argument("lane", nargs="+", help="file, name=file or name=1000,500")
    s.add_argument("-e", "--enzyme", action="append", help="digest with these enzymes first")
    s.add_argument("--out", default="gel.png")

    s = add("orf", "find ORFs")
    s.add_argument("file")
    s.add_argument("--min-aa", type=int, default=100)
    s.add_argument("--limit", type=int, default=20)

    s = add("translate", "translate to protein")
    s.add_argument("file")
    s.add_argument("start", type=int)
    s.add_argument("end", type=int)
    s.add_argument("--strand", type=int, default=1, choices=[1, -1])
    s.add_argument("--table", type=int, default=1)

    s = add("find", "search a sequence")
    s.add_argument("file")
    s.add_argument("query")
    s.add_argument("--limit", type=int, default=30)

    s = add("clone", "cloning")
    s.add_argument("file", nargs="+")
    s.add_argument("-m", "--method", required=True,
                   help="ligation|restriction|gibson|goldengate|infusion|gateway-bp|gateway-lr")
    s.add_argument("-e", "--enzyme", action="append")
    s.add_argument("--limit", type=int, default=25, help="Gibson minimum overlap")
    s.add_argument("--linear", action="store_true", help="treat inputs as linear")
    s.add_argument("--allow-linear", action="store_true", help="also allow linear products")
    s.add_argument("--out")

    s = add("pcr", "PCR product")
    s.add_argument("file")
    s.add_argument("fwd")
    s.add_argument("rev")
    s.add_argument("--limit", type=int, default=14, help="minimum primer match length")
    s.add_argument("--out")

    s = add("primers", "design primers")
    s.add_argument("file")
    s.add_argument("start", type=int)
    s.add_argument("end", type=int)
    s.add_argument("--tm", type=float, default=58.0)
    s.add_argument("--fwd-tail", default="", help="5' extension (e.g. a restriction site)")
    s.add_argument("--rev-tail", default="")

    s = add("check", "primer quality")
    s.add_argument("primer", nargs="+")

    s = add("fold", "secondary structure")
    s.add_argument("seq", help="sequence or file")

    s = add("annotate", "auto-detect features")
    s.add_argument("file")
    s.add_argument("--min-identity", type=float, default=annotate_mod.MIN_IDENTITY)
    s.add_argument("--keep", action="store_true", help="keep existing features")
    s.add_argument("--out")

    s = add("mutate", "mutation")
    s.add_argument("file")
    s.add_argument("action", choices=["sub", "ins", "del"])
    s.add_argument("start", type=int)
    s.add_argument("--end", type=int, help="end position for 'del'")
    s.add_argument("--new", help="new sequence for 'sub'/'ins'")
    s.add_argument("--primers", action="store_true", help="also output mutagenesis primers")
    s.add_argument("--out")

    s = add("origin", "move the origin")
    s.add_argument("file")
    s.add_argument("position", type=int)
    s.add_argument("--out")

    s = add("flip", "reverse complement")
    s.add_argument("file")
    s.add_argument("--out")

    s = add("feature", "add a feature")
    s.add_argument("file")
    s.add_argument("start", type=int)
    s.add_argument("end", type=int)
    s.add_argument("label")
    s.add_argument("--type", default="misc_feature")
    s.add_argument("--strand", type=int, default=1, choices=[1, -1])
    s.add_argument("--out")

    s = add("verify", "verify a sequencing read")
    s.add_argument("file", help="reference plasmid")
    s.add_argument("read", help="read (.ab1 or fasta/gb)")

    s = add("methyl", "methylation blocking")
    s.add_argument("file")
    s.add_argument("-e", "--enzyme", action="append")
    s.add_argument("--max-cuts", type=int)
    s.add_argument("--all-enzymes", action="store_true")

    s = add("guides", "CRISPR guides")
    s.add_argument("file")
    s.add_argument("--pam", default="NGG")
    s.add_argument("--length", type=int, default=20)
    s.add_argument("--gc-min", type=float)
    s.add_argument("--gc-max", type=float)
    s.add_argument("--limit", type=int, default=20)

    s = add("align", "multiple alignment")
    s.add_argument("file", nargs="+")
    s.add_argument("--program", default="clustalo", choices=["clustalo", "muscle"])
    s.add_argument("--out")

    add("doctor", "check the installation and list what works", aliases=("conf",))
    add("status", "what is available on this machine")
    return p


def _use_utf8_output() -> None:
    """Make stdout/stderr able to carry the characters the reports use.

    Windows defaults to the legacy code page (cp1252), which cannot encode
    the arrows, degree signs and check marks in the output, so printing
    raised UnicodeEncodeError and `sg conf` crashed there. ``errors=replace``
    keeps the process alive on any console that still cannot render a glyph.
    """
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, OSError, ValueError):
            pass  # redirected or already-wrapped stream: leave it alone


def main(argv=None) -> int:
    _use_utf8_output()
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.cmd:
        print(HELP)
        return 0
    call, render = COMMANDS[args.cmd]
    try:
        result = call(args)
    except SgError as e:
        if args.json:
            print(json.dumps({"error": str(e)}))
        else:
            print(f"Error: {e}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        return 130
    print(json.dumps(result, ensure_ascii=False, indent=2) if args.json
          else render(result, args))
    return 0 if args.cmd != "doctor" or result["ok"] else 2


if __name__ == "__main__":
    raise SystemExit(main())

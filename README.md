# sgserver

Command-line toolkit **and MCP server** for the everyday work you would do in
SnapGene: inspect plasmids, find restriction sites, digest, clone, run PCR,
design primers, edit sequences, auto-annotate features, verify sequencing
reads and find CRISPR guides.

It is built so that **humans and AI agents use the same tools**: every command
has `--json` output, and every operation is also exposed as an MCP tool.

> **Independent project.** SnapGene is a trademark of GSL Biotech LLC. This
> repository contains no SnapGene software, data or licence. Features that use
> SnapGene's own programs call *your* installation (see [SnapGene is
> optional](#snapgene-is-optional)).

## Install

```bash
# CLI (+ MCP server)
uv tool install "sgserver[mcp] @ git+https://github.com/ikoshos-gland/sgserver"
# or from a clone
git clone https://github.com/ikoshos-gland/sgserver && cd sgserver
uv tool install --editable ".[mcp]"
```

Requires Python ≥ 3.12. Then:

```bash
sg --help
sg conf                # check installation & active capabilities (alias: sg doctor)
sg --json conf
```

## Quick start

```bash
sg info plasmid.gb                          # length, topology, GC, features
sg sites plasmid.gb --max-cuts 1            # unique cutters
sg digest plasmid.gb -e BamHI -e EcoRV      # fragments and ends
sg gel "cut=plasmid.gb" -e BamHI --out gel.png
sg orf plasmid.gb --min-aa 100
sg translate plasmid.gb 1653 2372
sg find plasmid.gb GGATCC                   # IUPAC letters work
```

### Insert a gene into a vector (full example)

Directional cloning of mCherry into pUC19 with EcoRI and BamHI:

```bash
# 1) primers with restriction-site tails (GCGC helps the enzyme bind)
sg primers mCherry.gb 1 711 --fwd-tail GCGCGAATTC --rev-tail GCGCGGATCC
# 2) amplify → 731 bp with the sites at the ends
sg pcr mCherry.gb GCGCGAATTCATGGTGAGCAAGGGCGAG GCGCGGATCCCTACTTGTACAGCTCGTCC --out amp.gb
# 3) clone
sg clone pUC19.gb amp.gb -m restriction -e EcoRI -e BamHI --out pUC19-mCherry.gb
```

Result: 3,382 bp = 2,686 (pUC19) − 21 (removed fragment) + 717 (cut insert).
Both sites are recreated and the downstream SalI/HindIII sites survive.

> **Pitfall.** Do not write digested fragments to files and join them later:
> GenBank cannot store sticky ends, so the product comes out a few bp too long.
> Pass the whole plasmids and name the enzyme with `-e`.

## Commands

| Group | Commands |
|---|---|
| Read / convert | `info`, `convert`, `map`, `history` |
| Analysis | `sites`, `digest`, `gel`, `orf`, `translate`, `find` |
| Cloning / PCR | `clone`, `pcr`, `primers`, `check`, `fold` |
| Editing | `annotate`, `mutate`, `origin`, `flip`, `feature` |
| Validation | `verify`, `methyl`, `guides` |
| Alignment | `align` |
| Setup | `doctor` (alias: `conf`), `status` |

Cloning methods: `restriction`, `goldengate`, `gibson`, `infusion`,
`ligation`, `gateway-bp`, `gateway-lr`.

Positions are 1-based and inclusive; `strand` is `1` or `-1`.

## For AI agents

### JSON output

Put `--json` **before** the command:

```bash
sg --json digest plasmid.gb -e EcoRI -e BamHI
```
```json
{
  "enzymes": ["EcoRI", "BamHI"],
  "fragments": [
    {"index": 1, "size_bp": 2665, "end_5": "5' overhang GATC", "end_3": "5' overhang AATT"},
    {"index": 2, "size_bp": 21,   "end_5": "5' overhang AATT", "end_3": "5' overhang GATC"}
  ],
  "files": []
}
```

Errors are `{"error": "..."}` with exit code 1. The same dictionaries come from
the Python API in [`sgserver.api`](src/sgserver/api.py).

### MCP server

`sg-mcp` speaks MCP over stdio and exposes 25 tools (`sg_info`, `sg_sites`,
`sg_digest`, `sg_clone`, `sg_pcr`, `sg_design_primers`, `sg_mutate`,
`sg_verify_read`, ...). Add it to your client:

```json
{
  "mcpServers": {
    "sgserver": {
      "command": "sg-mcp",
      "env": { "SG_MCP_ROOT": "/path/to/your/project" }
    }
  }
}
```

- **Claude Code:** `claude mcp add sgserver -e SG_MCP_ROOT=$PWD -- sg-mcp`
- **Codex / Gemini / others:** use the JSON block above in their MCP settings.

**Safety:** every path an agent passes must resolve inside `SG_MCP_ROOT`
(default: the directory the server starts in); anything else is rejected. The
server uses local stdio only and opens no network port.

### Skill

[`skills/sgserver/SKILL.md`](skills/sgserver/SKILL.md) teaches an agent the
commands and the pitfalls. Copy the folder into your agent's skills directory
(e.g. `~/.claude/skills/`).

## SnapGene is optional

| Works without SnapGene | Needs a SnapGene installation |
|---|---|
| reading `.dna`, `.gb`, `.fasta`, `.embl`; writing `.gb`/`.fasta` | writing `.dna`, `convert` to/from `.dna` |
| `info`, `sites`, `digest`, `gel`, `orf`, `translate`, `find` | `map` (map image), `history` (SnapGene's own export) |
| `clone`, `pcr`, `primers`, `check`, `mutate`, `origin`, `flip`, `feature` | `annotate` (uses SnapGene's feature database) |
| `verify`, `methyl`, `guides` | `fold` (ViennaRNA), `align` (Clustal Omega / MUSCLE) |

SnapGene is looked up in `/opt/gslbiotech/snapgene`; set `SNAPGENE_DIR` to use
another location. The CLI-based features also need a valid SnapGene licence
(`~/.snapgene/snapgene-license.lic`). **Nothing from SnapGene is copied into
this repository.**

## Validation

The tests use SnapGene's own sample project as an independent reference: it
holds both the inputs of a cloning and the product SnapGene produced.

| Test | Result |
|---|---|
| pUC57 + mTAGBFP2, SalI/KpnI | 3,387 bp, sequence identical to SnapGene's |
| pGGA + mCherry, BsaI (Golden Gate) | 2,848 bp, sequence identical to SnapGene's |
| pUC19 feature detection | all 10 features SnapGene finds, no false positives |
| BamHI + EcoRV digest | 496 + 6,739 bp (fragments sum to plasmid size) |
| mEGFP translation | `MVSKGEE…` |

```bash
uv run python -m pytest tests -q
```

Without SnapGene installed, the tests that need its sample files are skipped
and the synthetic-sequence tests (`test_api.py`, `test_mcp.py`) still run — this
is what CI does.

## Known limits

- With the SnapGene window open, `convert`, `map` and `history` may silently
  do nothing; the error message reminds you.
- `history` uses SnapGene's own export and falls back to pydna's parser, which
  can choke on ambiguous letters (N, R, Y).
- Output files are GenBank. They open in SnapGene with their features, but are
  not recorded in SnapGene's own cloning history.
- Methylation: an overlap with Dam/Dcm is only a warning unless the enzyme's
  sensitivity is known (e.g. BamHI is unaffected by Dam, MboI is blocked, DpnI requires it).
- `.scf` trace files cannot be read (Biopython does not support them); export
  them as FASTA from SnapGene first.

## License

[MIT](LICENSE)

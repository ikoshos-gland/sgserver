---
name: sgserver
description: Work with DNA/plasmid sequence files like SnapGene does — inspect, find restriction sites, digest, clone (restriction, Golden Gate, Gibson), run in-silico PCR, design primers, edit sequences, auto-annotate features, verify sequencing reads, find CRISPR guides. Use when the user mentions plasmids, cloning, primers, digests, SnapGene, .dna/.gb/.fasta files, or "/conf".
---

# sgserver

## Step 0 — Run configuration check (`sg conf`)
Before your first command in a session (or whenever the user asks `/conf` or to check setup), verify what is installed and which capabilities are active on this machine:

- **MCP:** call `sg_doctor()`
- **CLI:**
  ```bash
  sg --json conf
  ```

`capabilities` in the output tells you whether SnapGene-specific features (`.dna` writing, `map`, `history`, `annotate`, `fold`, `align`) are available or if you should stick to the SnapGene-free core (`.gb`/`.fasta`, `info`, `sites`, `digest`, `clone`, `pcr`, `primers`, `mutate`, `verify`, `guides`, `methyl`, `gel`). If anything is missing, `next_steps` gives the exact fix command.

## Calling conventions
Prefer MCP tools if the `sgserver` MCP server is connected (`sg_doctor`, `sg_info`, `sg_digest`, `sg_clone`, ...); otherwise call the CLI with `--json` so results are machine-readable:

```bash
sg --json <command> <args>     # --json goes BEFORE the command
```

Errors come back as `{"error": "..."}` with exit code 1.

- Positions are **1-based and inclusive**. `strand` is `1` or `-1`.

## Common tasks
| Goal | CLI |
|---|---|
| Check setup & capabilities | `sg --json conf` (or `sg doctor`) |
| Overview | `sg --json info plasmid.gb` |
| Unique cutters | `sg --json sites plasmid.gb --max-cuts 1` |
| Digest | `sg --json digest plasmid.gb -e EcoRI -e BamHI` |
| Insert a gene (directional) | `sg --json primers gene.gb 1 711 --fwd-tail GCGCGAATTC --rev-tail GCGCGGATCC`, then `sg --json pcr gene.gb <fwd> <rev> --out amp.gb`, then `sg --json clone vector.gb amp.gb -m restriction -e EcoRI -e BamHI --out new.gb` |
| Golden Gate | `sg --json clone vec.gb part.gb -m goldengate -e BsaI --out new.gb` |
| Gibson | `sg --json clone a.gb b.gb c.gb -m gibson --out new.gb` |
| Point mutation + primers | `sg --json mutate p.gb sub 500 --new GG --primers --out m.gb` |
| Check a sequencing read | `sg --json verify reference.gb read.ab1` |
| CRISPR guides | `sg --json guides p.gb --gc-min 40 --gc-max 70` |
| Methylation overlap | `sg --json methyl p.gb` |

## Pitfalls — read before cloning
1. **Pass whole plasmids plus enzyme names to `clone`.** Never digest, write the
   fragments to files, and join them: GenBank cannot store sticky ends, so the
   product comes out a few bp too long.
2. `clone` can return several products; small ones are usually by-products.
   The wanted construct is normally the **largest**.
3. A methylation `status` of `overlap` is only a warning; `blocked` /
   `unaffected` / `required` are known cases (e.g. BamHI is *unaffected* by Dam).
4. If the SnapGene window is open, `convert`, `map` and `history` may silently
   do nothing — `sg conf` warns about this; ask the user to close it.
5. Auto-annotation uses SnapGene's feature database, so it needs SnapGene installed.

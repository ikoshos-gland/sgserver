# Contributing

```bash
git clone https://github.com/ikoshos-gland/sgserver && cd sgserver
uv sync --all-extras
uv run python -m pytest tests -q
```

- Add the logic to a module, expose it in `src/sgserver/api.py` (returns a
  JSON-serialisable `dict`), then wire the CLI command (`cli.py`) and the MCP
  tool (`mcp_server.py`). Keep the three in sync.
- Positions are 1-based and inclusive. Raise `SgError` for expected errors.
- Tests that need SnapGene's sample files use `requires_samples`; prefer tests
  built on synthetic sequences (see `tests/test_api.py`) so CI can run them.
- **Never commit** SnapGene files, licence files (`*.lic`) or sample data.
EOF

# Contributing to RAS Commander Hydro

Thank you for your interest in contributing to **ras-commander-hydro** -- Python tools for hydrologic analysis workflows built on ArcHydro patterns. This project is maintained by [CLB Engineering Corporation](https://clbengineering.com/) and released under the MIT license.

## Our Philosophy: Don't Ask Me, Ask a GPT!

This project was built with LLMs, and we welcome LLM-assisted contributions. Use whatever agent or model works for you -- Claude Code, Cursor, Copilot, Aider, ChatGPT, or anything else. There is no gatekeeping on tooling.

The tradeoff is simple: **if an LLM helped you write it, an LLM should also help you review it.** We call this the Self-Review Contract. LLM self-review catches the majority of issues before a maintainer ever looks at your PR, which keeps the review cycle fast for everyone.

This approach follows the [LLM Forward](https://clbengineering.com/llm-forward) philosophy pioneered by CLB Engineering -- technology accelerates engineering insight without replacing professional judgment, especially in safety-critical domains like water resources.

---

## Quick Start

### 1. Fork and Clone

```bash
git clone https://github.com/<your-username>/ras-commander-hydro.git
cd ras-commander-hydro
```

### 2. Install Dependencies

```bash
pip install -e .
```

Or install dependencies manually based on your workflow. Core dependencies include GDAL/OGR, geopandas, rasterio, shapely, and numpy.

### 3. Create a Branch

```bash
git checkout -b feat/your-feature-name
```

### 4. Launch Your Agent

Open the repo in your preferred LLM-assisted development environment. Point it at this file and the README for context. Start building.

### 5. Self-Review and Submit

Before opening a PR, run through the LLM Self-Review Checklist below. Then push your branch and open a pull request using the provided template.

---

## The Self-Review Contract

Every pull request should include evidence that the contributor (or their LLM) reviewed the changes against the checklist below. This is not bureaucracy -- it is a practical way to reduce back-and-forth during review.

In your PR description, include the LLM Self-Review section from the pull request template. Check off the items that apply. If something does not apply, note why.

Maintainers trust that you ran the checklist honestly. In return, reviews focus on design and domain questions rather than style nits.

---

## LLM Self-Review Checklist

Ask your LLM to review your changes against these criteria before you open the PR.

### Code Quality

- [ ] All public functions have docstrings with Args, Returns, and Raises sections
- [ ] Functions use `@log_call` or appropriate logging for traceability
- [ ] Error handling is explicit -- no bare `except:` clauses
- [ ] File paths use `pathlib.Path`, not string concatenation
- [ ] No hardcoded absolute paths -- accept paths as parameters
- [ ] Type hints on function signatures where practical

### Hydrology Specifics

- [ ] ArcHydro workflows follow established naming conventions (e.g., FDR, FAC, STR)
- [ ] Geospatial data handling preserves CRS throughout the pipeline
- [ ] CRS transformations are explicit, not silent
- [ ] Raster operations handle NoData values correctly
- [ ] Vector operations validate geometry before processing
- [ ] Units are documented in docstrings (feet, meters, cfs, cms)

### Testing

- [ ] Tested with real geospatial data, not synthetic mocks
- [ ] Edge cases considered (empty datasets, single-feature inputs, large rasters)
- [ ] Results are visually verifiable where appropriate (maps, plots, tables)
- [ ] No test relies on network access unless explicitly marked

### Documentation

- [ ] New functions are documented in relevant module docstrings
- [ ] Example usage included in docstring or notebook
- [ ] Breaking changes noted in the PR description

---

## What We Accept

- **Hydrology tools**: Watershed delineation, flow accumulation, stream network extraction, catchment analysis, and related ArcHydro-pattern workflows
- **ArcHydro extensions**: New processing steps that fit the ArcHydro pipeline
- **Geospatial utilities**: Raster/vector helpers for hydrologic analysis (DEM processing, catchment polygons, pour point snapping)
- **Integration with ras-commander**: Tools that bridge hydrologic analysis with HEC-RAS modeling
- **Documentation improvements**: Clearer docstrings, better examples, notebook tutorials
- **Bug fixes**: Especially those tested against real data

## What We Don't Accept

- **Breaking API changes without discussion**: Open an issue first to discuss the motivation and migration path
- **Unjustified new dependencies**: Every dependency adds maintenance burden. If you need a new package, explain why existing tools cannot accomplish the task
- **Mock-based tests as the sole testing strategy**: We test with real geospatial data. Mocks can supplement but should not replace real-data tests
- **Code without docstrings**: Public functions must be documented

---

## Commit Messages

We use [Conventional Commits](https://www.conventionalcommits.org/):

```
feat(watershed): Add pour point snapping to stream network
fix(fac): Handle NoData cells in flow accumulation grid
docs(examples): Add catchment delineation notebook
refactor(dem): Simplify pit-filling pipeline
test(stream): Add real-data test for stream order extraction
```

### LLM Attribution

If an LLM assisted with your contribution, include a `Co-Authored-By` trailer:

```
feat(catchment): Add upstream area calculation

Implement recursive upstream area accumulation from flow direction grid.

Co-Authored-By: Claude <noreply@anthropic.com>
```

This is encouraged, not required. Attribution helps the community understand how LLMs are being used in practice.

---

## Branch Naming

Use descriptive prefixes:

- `feat/` -- New features
- `fix/` -- Bug fixes
- `docs/` -- Documentation only
- `refactor/` -- Code restructuring without behavior change
- `test/` -- Adding or improving tests

---

## Community Standards

### Safety-Critical Domain

Hydrologic analysis supports flood risk management, infrastructure design, and water resources planning. These are safety-critical applications where errors can affect public welfare. Contributions should reflect this responsibility:

- Validate results against known benchmarks where possible
- Document assumptions and limitations clearly
- Prefer conservative defaults over aggressive optimizations
- When in doubt, raise an error rather than return a plausible-but-wrong result

### Professional Conduct

Be respectful, constructive, and patient. We welcome contributors of all experience levels. If someone's PR needs work, explain what and why -- don't just reject it.

### LLM Forward

We believe LLMs are powerful tools for engineering when used with professional judgment and proper verification. Read more about the LLM Forward approach at [clbengineering.com/llm-forward](https://clbengineering.com/llm-forward).

---

## Getting Help

- **Open an issue**: For bugs, feature requests, or questions about the codebase
- **Start a discussion**: For broader topics about hydrologic workflows or ArcHydro patterns
- **Read the README**: For an overview of the project and its capabilities
- **Ask your LLM**: Seriously -- point it at this repo and ask it to explain the code. That is what it is for.

---

*Maintained by [CLB Engineering Corporation](https://clbengineering.com/). MIT License.*

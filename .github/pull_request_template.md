## Summary

<!-- Describe what this PR does in 2-3 sentences. -->

## Type of Change

- [ ] New feature (non-breaking addition)
- [ ] Bug fix (non-breaking fix for an existing issue)
- [ ] Documentation (docstrings, notebooks, or markdown)
- [ ] Refactor (no behavior change)
- [ ] Breaking change (requires migration -- explain below)

## LLM Self-Review

<!-- Check off items you or your LLM verified. Mark N/A if not applicable. -->

### Code Quality
- [ ] Public functions have docstrings (Args, Returns, Raises)
- [ ] Logging and error handling are explicit (no bare `except:`)
- [ ] File paths use `pathlib.Path`
- [ ] Type hints on function signatures

### Hydrology Specifics
- [ ] CRS is preserved and transformations are explicit
- [ ] NoData and edge cases are handled
- [ ] Units documented in docstrings

### Testing
- [ ] Tested with real geospatial data
- [ ] Edge cases considered
- [ ] Results are visually verifiable where appropriate

## Test Plan

<!-- Describe how you tested these changes. Include data sources, commands, or notebooks used. -->

- [ ] Ran existing tests (`python run_tests.py`)
- [ ] Tested manually with real data (describe below)
- [ ] Added new tests

## LLM Attribution

<!-- If an LLM assisted, note which one. This is encouraged, not required. -->

- [ ] LLM-assisted (model: _______________)
- [ ] No LLM assistance

## Additional Notes

<!-- Any context, screenshots, or references that help reviewers. -->

# Testing Guide — RAS Commander Arc Hydro Tools

## Quick Start

Run tests using ArcGIS Pro's Python environment:

```batch
run_tests.bat tier1          # Local testdata/ — fast, always available
run_tests.bat tier2          # Production models on L:\ drive
run_tests.bat all            # Everything
run_tests.bat unicode        # Unicode/encoding tests only
run_tests.bat unit           # Unit tests (mostly no arcpy needed)
run_tests.bat integration    # Integration tests (requires arcpy)
```

Or directly via Python:

```
python run_tests.py tier1
python run_tests.py tier1 -- -x --pdb    # Stop on first failure, drop into debugger
```

## Test Tiers

### Tier 1 — Local Test Data (`testdata/`)

Always available. Uses 4 HDF files covering 1D steady, 1D unsteady, 2D unsteady, and 2D with pipe networks.

| Test File | HDF Source | What It Tests |
|-----------|-----------|---------------|
| BaldEagle.p01.hdf | 1D Unsteady | Cross sections, centerlines, bank/edge lines |
| BaldEagle.p02.hdf | 1D Steady | Cross sections without unsteady results |
| BaldEagleDamBrk.p07.hdf | 2D Unsteady | Mesh geometry, BC lines, breaklines, Max WSE/velocity |
| DavisStormSystem.p02.hdf | 2D + Pipes | Pipe conduits, pipe nodes, pipe networks |

### Tier 2 — Production Models (`L:\`)

Auto-discovers `.p*.hdf` files across `L:\Region_1` through `L:\Region_7`. Automatically skipped if the `L:\` drive is not mounted.

Tier 2 tests:
- **Strict UTF-8 scan**: Decodes ALL string attributes with strict UTF-8 (no `'ignore'`) to surface encoding issues that production code silently handles
- **Structure validation**: Verifies every file opens and has expected groups
- **Tool execution**: Runs `OrganizeRASProject` on a random sample of 5 files per region

## Test Categories

### Unit Tests (`tests/unit/`)

No arcpy required for most tests. Uses only h5py and numpy.

- **test_utils.py** — Pure function tests: `extract_project_and_plan_info`, `get_feature_dataset_name`, `is_clockwise_numpy`, `get_dynamic_fields_from_data`
- **test_hdf_reading.py** — HDF structure validation: expected groups, datasets, attributes
- **test_unicode_handling.py** — Byte decode edge cases, verifies bug fixes in `rc_utils.py`

### Integration Tests (`tests/integration/`)

Require arcpy. Execute tools end-to-end using `MockParam`/`MockMessages` pattern.

- **test_load_1d_geometry.py** — Cross sections, centerlines, bank lines, structures
- **test_load_2d_geometry.py** — Mesh perimeters, cell centers, faces, breaklines, BC lines, pipe networks
- **test_load_2d_results.py** — Max WSE, max velocity, error handling for 1D files
- **test_organize_ras_project.py** — Batch processing, directory scanning, `load_to_map=False`

### Edge Case Tests (`tests/edge_cases/`)

- **test_missing_hdf_groups.py** — Wrong tool on wrong file type, empty HDF files, `_check_hdf_contents` validation
- **test_unicode_paths.py** — Non-ASCII directory names, unicode filenames, spaces in paths
- **test_empty_datasets.py** — Zero-length HDF arrays, all-NaN attribute tables, missing expected datasets

## Pytest Markers

| Marker | Description |
|--------|-------------|
| `tier1` | Uses local testdata/ (always available) |
| `tier2` | Uses production models on L:\ (auto-skip if unavailable) |
| `requires_arcpy` | Needs ArcGIS Pro Python environment |
| `slow` | Takes more than 30 seconds |
| `unicode` | Unicode/encoding-specific tests |

Run specific markers:

```batch
run_tests.bat tier1 -m unicode              # Unicode tests from tier1 only
run_tests.bat tier1 -m "not requires_arcpy" # Unit tests without arcpy
run_tests.bat tier1 -m "not slow"           # Skip slow tests
```

## MockParam / MockMessages Pattern

Tests use lightweight stand-ins for arcpy objects defined in `tests/conftest.py`:

```python
from conftest import MockParam, MockMessages

# Build parameter list matching tool's getParameterInfo() indices
params = [
    MockParam(hdf_path, hdf_path),  # value and valueAsText
    MockParam(None),                 # optional param
    MockParam(["Cross Sections"]),   # multi-value param (sets .values)
    ...
]

messages = MockMessages()
tool.execute(params, messages)

assert not messages.has_error()
```

## Bug Fixes Verified by Tests

Two strict UTF-8 decode bugs were fixed in `Scripts/archydro/rc_utils.py`:

1. **Line 32**: `proj_wkt_attr.decode("utf-8")` → `proj_wkt_attr.decode("utf-8", "ignore")`
2. **Line 290**: `time_str.decode('utf-8')` → `time_str.decode('utf-8', 'ignore')`

The `test_unicode_handling.py::TestBugFixVerification` class inspects the source code to confirm these fixes remain in place.

## Adding New Tests

1. Place test files in the appropriate category directory
2. Use `pytestmark = [pytest.mark.tier1]` (or tier2) at module level
3. Use fixtures from `tests/conftest.py` for HDF file paths and mock objects
4. For integration tests, build parameter lists matching the tool's `getParameterInfo()` index layout
5. Output to `memory\` workspace for in-memory feature classes (no cleanup needed)

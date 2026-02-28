# -*- coding: utf-8 -*-
"""
Tier 2 conftest - Production model discovery on L:\\ drive.

Auto-discovers HDF files across L:\\Region_1 through L:\\Region_7.
Skips all tier2 tests if L:\\ is not mounted/accessible.
"""

import os
import pytest
import random

L_DRIVE_BASE = r"L:\\"
REGION_DIRS = [f"Region_{i}" for i in range(1, 8)]
MAX_DEPTH = 4
HDF_PATTERN_SUFFIXES = (".hdf",)


def _is_plan_hdf(filename):
    """Check if filename matches the p*.hdf pattern."""
    import re
    return bool(re.search(r"\.p\d{2}\.hdf$", filename, re.IGNORECASE))


def discover_hdf_files(base_path=L_DRIVE_BASE, max_depth=MAX_DEPTH):
    """Walk L:\\ regions and collect all .p*.hdf files.

    Returns:
        dict: {region_name: [list of hdf paths]}
    """
    results = {}

    for region_dir in REGION_DIRS:
        region_path = os.path.join(base_path, region_dir)
        if not os.path.isdir(region_path):
            continue

        hdf_files = []
        for root, dirs, files in os.walk(region_path):
            # Limit depth
            depth = root.replace(region_path, "").count(os.sep)
            if depth >= max_depth:
                dirs.clear()
                continue

            for f in files:
                if _is_plan_hdf(f):
                    hdf_files.append(os.path.join(root, f))

        if hdf_files:
            results[region_dir] = sorted(hdf_files)

    return results


def _l_drive_available():
    """Check if L:\\ drive is mounted and has at least one region dir."""
    if not os.path.isdir(L_DRIVE_BASE):
        return False
    for region_dir in REGION_DIRS:
        if os.path.isdir(os.path.join(L_DRIVE_BASE, region_dir)):
            return True
    return False


# ---------------------------------------------------------------------------
# Session-scoped discovery
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session")
def production_hdf_map():
    """Discover all production HDF files, skip if L:\\ unavailable."""
    if not _l_drive_available():
        pytest.skip("L:\\ drive not available — skipping tier2 tests")
    hdf_map = discover_hdf_files()
    if not hdf_map:
        pytest.skip("No HDF files found on L:\\ drive")
    return hdf_map


@pytest.fixture(scope="session")
def all_production_hdfs(production_hdf_map):
    """Flat list of all discovered production HDF paths."""
    all_files = []
    for files in production_hdf_map.values():
        all_files.extend(files)
    return all_files


@pytest.fixture(scope="session")
def sampled_production_hdfs(production_hdf_map):
    """Random sample of up to 5 HDF files per region for tool execution tests."""
    sampled = []
    for region, files in production_hdf_map.items():
        sample_size = min(5, len(files))
        sampled.extend(random.sample(files, sample_size))
    return sampled

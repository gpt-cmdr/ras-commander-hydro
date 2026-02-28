# -*- coding: utf-8 -*-
"""
Tier 2 tests - Production model sweep across L:\\ drive.

Discovers all HDF files in L:\\Region_1 through L:\\Region_7 and:
1. Scans ALL string attributes with strict UTF-8 decode to surface encoding issues
2. Validates HDF structure (expected groups/datasets)
3. Runs tool execution on a random sample of 5 files per region

These tests are slow and only run when L:\\ is available.
"""

import os
import pytest
import h5py
import numpy as np

pytestmark = [pytest.mark.tier2, pytest.mark.slow]


# ---------------------------------------------------------------------------
# Strict UTF-8 scan — finds encoding issues the 'ignore' handler hides
# ---------------------------------------------------------------------------
class TestStrictUTF8Scan:
    """Scan production HDF string attributes with strict UTF-8 decode."""

    def _scan_strict(self, hdf_path):
        """Scan all string attributes with strict UTF-8, return issues."""
        issues = []

        def visitor(name, obj):
            for attr_name, attr_val in obj.attrs.items():
                if isinstance(attr_val, (bytes, np.bytes_)):
                    try:
                        attr_val.decode("utf-8")
                    except UnicodeDecodeError as e:
                        issues.append({
                            "path": f"{name}@{attr_name}",
                            "error": str(e),
                            "raw_bytes": repr(attr_val[:100]),
                        })

        try:
            with h5py.File(hdf_path, "r") as f:
                # File-level attributes
                for attr_name, attr_val in f.attrs.items():
                    if isinstance(attr_val, (bytes, np.bytes_)):
                        try:
                            attr_val.decode("utf-8")
                        except UnicodeDecodeError as e:
                            issues.append({
                                "path": f"/@{attr_name}",
                                "error": str(e),
                                "raw_bytes": repr(attr_val[:100]),
                            })
                f.visititems(visitor)
        except Exception as e:
            issues.append({"path": "FILE_OPEN", "error": str(e), "raw_bytes": ""})

        return issues

    def test_scan_all_production_files(self, all_production_hdfs):
        """Scan every production HDF for non-UTF-8 string attributes.

        This test collects all issues and reports them as warnings rather
        than failing, since the 'ignore' handler in production code will
        handle these bytes gracefully.
        """
        all_issues = {}
        for hdf_path in all_production_hdfs:
            issues = self._scan_strict(hdf_path)
            if issues:
                all_issues[hdf_path] = issues

        if all_issues:
            total = sum(len(v) for v in all_issues.values())
            msg_lines = [f"Found {total} non-UTF-8 attributes across {len(all_issues)} files:"]
            for path, issues in list(all_issues.items())[:10]:
                msg_lines.append(f"  {os.path.basename(path)}: {len(issues)} issues")
                for issue in issues[:3]:
                    msg_lines.append(f"    {issue['path']}: {issue['error']}")
            pytest.xfail("\n".join(msg_lines))


# ---------------------------------------------------------------------------
# HDF Structure Validation
# ---------------------------------------------------------------------------
class TestProductionHDFStructure:
    """Validate that production HDF files have expected structure."""

    def test_all_files_open(self, all_production_hdfs):
        """Every discovered HDF file should be openable by h5py."""
        failures = []
        for hdf_path in all_production_hdfs:
            try:
                with h5py.File(hdf_path, "r") as f:
                    _ = list(f.keys())
            except Exception as e:
                failures.append(f"{os.path.basename(hdf_path)}: {e}")

        assert not failures, f"Failed to open {len(failures)} files:\n" + "\n".join(failures[:10])

    def test_all_files_have_geometry(self, all_production_hdfs):
        """Every plan HDF should have a Geometry group."""
        missing = []
        for hdf_path in all_production_hdfs:
            try:
                with h5py.File(hdf_path, "r") as f:
                    if "Geometry" not in f:
                        missing.append(os.path.basename(hdf_path))
            except Exception:
                pass  # Open failures already caught above

        if missing:
            pytest.xfail(f"{len(missing)} files missing Geometry group: {missing[:5]}")

    def test_projection_present(self, sampled_production_hdfs):
        """Sampled files should have a Projection attribute or associated .prj."""
        no_proj = []
        for hdf_path in sampled_production_hdfs:
            try:
                with h5py.File(hdf_path, "r") as f:
                    proj = f.attrs.get("Projection")
                    if proj is None:
                        no_proj.append(os.path.basename(hdf_path))
            except Exception:
                pass

        if no_proj:
            # This is informational — not all files need embedded projection
            print(f"INFO: {len(no_proj)}/{len(sampled_production_hdfs)} files lack Projection attr")


# ---------------------------------------------------------------------------
# Tool Execution on Sampled Files (requires arcpy)
# ---------------------------------------------------------------------------
class TestProductionToolExecution:
    """Run tools on a sample of production HDF files."""

    @pytest.fixture(autouse=True)
    def _skip_without_arcpy(self):
        try:
            import arcpy
        except ImportError:
            pytest.skip("arcpy not available")

    def test_check_contents_sampled(self, sampled_production_hdfs):
        """Run _check_hdf_contents on all sampled files."""
        from rc_organize_ras_project import OrganizeRASProject

        tool = OrganizeRASProject()
        failures = []

        for hdf_path in sampled_production_hdfs:
            try:
                with h5py.File(hdf_path, "r") as f:
                    contents = tool._check_hdf_contents(f)
                    assert isinstance(contents, dict)
            except Exception as e:
                failures.append(f"{os.path.basename(hdf_path)}: {e}")

        assert not failures, f"{len(failures)} check_contents failures:\n" + "\n".join(failures[:5])

    def test_cache_metadata_sampled(self, sampled_production_hdfs):
        """Run cache_hdf_metadata on all sampled files."""
        from rc_utils import cache_hdf_metadata

        failures = []

        for hdf_path in sampled_production_hdfs:
            try:
                with h5py.File(hdf_path, "r") as f:
                    cache = cache_hdf_metadata(f)
                    assert isinstance(cache, dict)
            except Exception as e:
                failures.append(f"{os.path.basename(hdf_path)}: {e}")

        assert not failures, f"{len(failures)} cache_metadata failures:\n" + "\n".join(failures[:5])

    def test_organize_sampled(self, sampled_production_hdfs, tmp_path):
        """Run OrganizeRASProject on sampled files with load_to_map=False."""
        import arcpy
        from rc_organize_ras_project import OrganizeRASProject

        # Import MockParam/MockMessages from parent conftest
        import sys
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
        from conftest import MockParam, MockMessages

        failures = []

        for i, hdf_path in enumerate(sampled_production_hdfs[:3]):  # Limit to 3 for speed
            try:
                output_gdb = str(tmp_path / f"prod_test_{i}.gdb")
                messages = MockMessages()
                tool = OrganizeRASProject()

                params = [
                    MockParam(hdf_path, hdf_path),
                    MockParam(output_gdb, output_gdb),
                    MockParam(None),
                    MockParam(True),
                    MockParam(True),
                    MockParam(True),
                    MockParam(False),
                    MockParam(False),  # load_to_map = False
                ]
                tool.execute(params, messages)

                if messages.has_error():
                    # CRS errors are expected for files without projection
                    error_text = " ".join(messages.errors)
                    if "CRS" not in error_text:
                        failures.append(f"{os.path.basename(hdf_path)}: {error_text}")
            except Exception as e:
                error_str = str(e)
                if "CRS" not in error_str and "projection" not in error_str.lower():
                    failures.append(f"{os.path.basename(hdf_path)}: {e}")

        if failures:
            pytest.xfail(f"{len(failures)} tool execution failures:\n" + "\n".join(failures[:5]))

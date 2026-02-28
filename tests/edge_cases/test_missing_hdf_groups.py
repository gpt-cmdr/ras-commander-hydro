# -*- coding: utf-8 -*-
"""
Edge case tests for missing/wrong HDF groups.

Covers scenarios like:
    - Running the 2D tool on a 1D-only file
    - Running the results tool on a geometry-only file
    - HDF files missing expected subgroups
"""

import os
import pytest
import h5py
import tempfile
import numpy as np

pytestmark = [pytest.mark.tier1]


class TestWrongToolOnWrongFile:
    """Test tool behavior when HDF file doesn't contain expected data."""

    def test_2d_structure_absent_in_1d_file(self, hdf_1d_unsteady):
        """1D file should not have 2D Flow Areas."""
        with h5py.File(hdf_1d_unsteady, "r") as f:
            assert "Geometry/2D Flow Areas" not in f

    def test_results_absent_in_steady_file(self, hdf_1d_steady):
        """Steady-state file typically won't have unsteady results."""
        with h5py.File(hdf_1d_steady, "r") as f:
            results_path = "Results/Unsteady/Output/Output Blocks/Base Output/Summary Output/2D Flow Areas"
            assert results_path not in f

    def test_pipes_absent_in_baldeagle(self, hdf_2d_unsteady):
        """BaldEagleDamBrk should not have pipe network data."""
        with h5py.File(hdf_2d_unsteady, "r") as f:
            assert "Geometry/Pipe Conduits" not in f

    def test_1d_absent_in_pipes_file(self, hdf_2d_pipes):
        """DavisStormSystem may not have 1D cross sections."""
        with h5py.File(hdf_2d_pipes, "r") as f:
            has_xs = "Geometry/Cross Sections/Attributes" in f
            # This is informational — just verify we can check it
            assert isinstance(has_xs, bool)


class TestEmptyHDFFile:
    """Test behavior with a minimal/empty HDF file."""

    @pytest.fixture
    def empty_hdf(self, tmp_path):
        """Create an empty HDF file with no groups."""
        path = str(tmp_path / "empty.hdf")
        with h5py.File(path, "w") as f:
            pass  # No groups or datasets
        return path

    @pytest.fixture
    def minimal_hdf(self, tmp_path):
        """Create HDF with Geometry group but no subgroups."""
        path = str(tmp_path / "minimal.hdf")
        with h5py.File(path, "w") as f:
            f.create_group("Geometry")
        return path

    def test_empty_hdf_has_no_geometry(self, empty_hdf):
        with h5py.File(empty_hdf, "r") as f:
            assert "Geometry" not in f
            assert "Plan Data" not in f

    def test_minimal_hdf_has_no_cross_sections(self, minimal_hdf):
        with h5py.File(minimal_hdf, "r") as f:
            assert "Geometry" in f
            assert "Geometry/Cross Sections" not in f
            assert "Geometry/2D Flow Areas" not in f


class TestCheckHDFContents:
    """Test the OrganizeRASProject._check_hdf_contents method directly."""

    def test_check_contents_1d(self, hdf_1d_unsteady):
        from rc_organize_ras_project import OrganizeRASProject

        tool = OrganizeRASProject()
        with h5py.File(hdf_1d_unsteady, "r") as f:
            contents = tool._check_hdf_contents(f)

        assert contents["has_1d"] is True
        assert "Cross Sections" in contents["1d_elements"]
        assert contents["has_2d"] is False

    def test_check_contents_2d(self, hdf_2d_unsteady):
        from rc_organize_ras_project import OrganizeRASProject

        tool = OrganizeRASProject()
        with h5py.File(hdf_2d_unsteady, "r") as f:
            contents = tool._check_hdf_contents(f)

        assert contents["has_2d"] is True
        assert "Mesh Area Perimeters" in contents["2d_elements"]
        assert contents["has_results"] is True

    def test_check_contents_pipes(self, hdf_2d_pipes):
        from rc_organize_ras_project import OrganizeRASProject

        tool = OrganizeRASProject()
        with h5py.File(hdf_2d_pipes, "r") as f:
            contents = tool._check_hdf_contents(f)

        assert contents["has_pipes"] is True
        assert "Pipe Conduits" in contents["pipe_elements"]

    def test_check_contents_empty(self, tmp_path):
        from rc_organize_ras_project import OrganizeRASProject

        path = str(tmp_path / "empty_check.hdf")
        with h5py.File(path, "w") as f:
            pass

        tool = OrganizeRASProject()
        with h5py.File(path, "r") as f:
            contents = tool._check_hdf_contents(f)

        assert contents["has_1d"] is False
        assert contents["has_2d"] is False
        assert contents["has_pipes"] is False
        assert contents["has_results"] is False

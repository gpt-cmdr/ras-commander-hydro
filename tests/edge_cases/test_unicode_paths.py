# -*- coding: utf-8 -*-
"""
Edge case tests for non-ASCII directory and file paths.

Tests that HDF file reading works when files are located in directories
with unicode characters (accented letters, CJK, etc.).
"""

import os
import shutil
import pytest
import h5py
import numpy as np

pytestmark = [pytest.mark.tier1, pytest.mark.unicode]


class TestUnicodeDirectoryPaths:
    """Test HDF reading from directories with non-ASCII names."""

    @pytest.fixture
    def unicode_dir_hdf(self, hdf_1d_unsteady, tmp_path):
        """Copy a test HDF file into a directory with unicode chars."""
        unicode_dir = tmp_path / "Prüfung_Straße_测试"
        unicode_dir.mkdir(parents=True, exist_ok=True)
        dest = unicode_dir / "BaldEagle.p01.hdf"
        shutil.copy2(hdf_1d_unsteady, str(dest))
        return str(dest)

    def test_h5py_opens_unicode_path(self, unicode_dir_hdf):
        """h5py should open files in unicode-named directories."""
        with h5py.File(unicode_dir_hdf, "r") as f:
            assert "Geometry" in f

    def test_read_attributes_from_unicode_path(self, unicode_dir_hdf):
        """Should read HDF attributes from unicode path without errors."""
        with h5py.File(unicode_dir_hdf, "r") as f:
            if "Geometry/Cross Sections/Attributes" in f:
                attrs = f["Geometry/Cross Sections/Attributes"][()]
                assert len(attrs) > 0

    def test_projection_from_unicode_path(self, unicode_dir_hdf):
        """get_ras_projection_wkt should handle unicode paths."""
        from rc_utils import get_ras_projection_wkt

        # This may return None if projection is in a .prj file not copied,
        # but should not crash
        try:
            result = get_ras_projection_wkt(unicode_dir_hdf)
        except UnicodeError:
            pytest.fail("UnicodeError when reading projection from unicode path")


class TestUnicodeFilenames:
    """Test with HDF files that have unicode characters in filename."""

    @pytest.fixture
    def unicode_filename_hdf(self, hdf_1d_unsteady, tmp_path):
        """Copy HDF file with a unicode filename."""
        dest = tmp_path / "Flüsse_河流.p01.hdf"
        shutil.copy2(hdf_1d_unsteady, str(dest))
        return str(dest)

    def test_h5py_opens_unicode_filename(self, unicode_filename_hdf):
        with h5py.File(unicode_filename_hdf, "r") as f:
            assert "Geometry" in f

    def test_extract_project_info_unicode(self, unicode_filename_hdf):
        """extract_project_and_plan_info should handle unicode filenames."""
        from rc_utils import extract_project_and_plan_info

        proj, plan, base = extract_project_and_plan_info(unicode_filename_hdf)
        assert plan == "01"
        assert "Flüsse_河流" in proj


class TestSpacesInPaths:
    """Test with spaces in directory names (common on Windows)."""

    @pytest.fixture
    def spaced_path_hdf(self, hdf_1d_unsteady, tmp_path):
        """Copy HDF to a directory with spaces."""
        spaced_dir = tmp_path / "My HEC-RAS Models" / "Project Files"
        spaced_dir.mkdir(parents=True, exist_ok=True)
        dest = spaced_dir / "BaldEagle.p01.hdf"
        shutil.copy2(hdf_1d_unsteady, str(dest))
        return str(dest)

    def test_h5py_opens_spaced_path(self, spaced_path_hdf):
        with h5py.File(spaced_path_hdf, "r") as f:
            assert "Geometry" in f

    def test_project_info_from_spaced_path(self, spaced_path_hdf):
        from rc_utils import extract_project_and_plan_info

        proj, plan, base = extract_project_and_plan_info(spaced_path_hdf)
        assert proj == "BaldEagle"
        assert plan == "01"

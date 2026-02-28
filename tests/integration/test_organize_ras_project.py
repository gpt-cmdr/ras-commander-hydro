# -*- coding: utf-8 -*-
"""
Integration tests for OrganizeRASProject tool (end-to-end batch processing).

Requires arcpy (ArcGIS Pro Python environment).

Parameter layout (8 params, including new load_to_map):
    0: input_path (folder or file)
    1: output_gdb
    2: override_crs
    3: include_1d_geometry (bool)
    4: include_2d_geometry (bool)
    5: include_2d_results (bool)
    6: include_cell_polygons (bool)
    7: load_to_map (bool) — NEW, set False for testing
"""

import os
import pytest

pytestmark = [pytest.mark.tier1, pytest.mark.requires_arcpy, pytest.mark.slow]

import arcpy
if getattr(arcpy, "_is_mock", False):
    pytest.skip("arcpy is mocked — real arcpy required", allow_module_level=True)

from conftest import MockParam, MockMessages


def _build_organize_params(input_path, output_gdb, include_1d=True, include_2d=True,
                           include_results=True, include_polys=False, load_to_map=False):
    """Build 8-element parameter list for OrganizeRASProject."""
    return [
        MockParam(input_path, input_path),  # 0: input_path
        MockParam(output_gdb, output_gdb),  # 1: output_gdb
        MockParam(None),                     # 2: override_crs
        MockParam(include_1d),               # 3: include_1d_geometry
        MockParam(include_2d),               # 4: include_2d_geometry
        MockParam(include_results),          # 5: include_2d_results
        MockParam(include_polys),            # 6: include_cell_polygons
        MockParam(load_to_map),              # 7: load_to_map (False for testing)
    ]


class TestOrganizeSingleHDF:
    """Test organizing a single HDF file."""

    def test_organize_2d_unsteady(self, hdf_2d_unsteady, mock_messages, tmp_path):
        from rc_organize_ras_project import OrganizeRASProject

        output_gdb = str(tmp_path / "Organized2D.gdb")
        tool = OrganizeRASProject()
        params = _build_organize_params(
            hdf_2d_unsteady, output_gdb,
            include_1d=False, include_2d=True, include_results=True,
        )
        tool.execute(params, mock_messages)

        assert not mock_messages.has_error(), f"Errors: {mock_messages.errors}"
        assert arcpy.Exists(output_gdb), "Output GDB not created"

    def test_organize_1d_unsteady(self, hdf_1d_unsteady, mock_messages, tmp_path):
        from rc_organize_ras_project import OrganizeRASProject

        output_gdb = str(tmp_path / "Organized1D.gdb")
        tool = OrganizeRASProject()
        params = _build_organize_params(
            hdf_1d_unsteady, output_gdb,
            include_1d=True, include_2d=False, include_results=False,
        )
        tool.execute(params, mock_messages)

        assert not mock_messages.has_error(), f"Errors: {mock_messages.errors}"

    def test_organize_with_pipes(self, hdf_2d_pipes, mock_messages, tmp_path):
        from rc_organize_ras_project import OrganizeRASProject

        output_gdb = str(tmp_path / "OrganizedPipes.gdb")
        tool = OrganizeRASProject()
        params = _build_organize_params(
            hdf_2d_pipes, output_gdb,
            include_1d=False, include_2d=True, include_results=True,
        )
        tool.execute(params, mock_messages)

        assert not mock_messages.has_error(), f"Errors: {mock_messages.errors}"


class TestOrganizeDirectory:
    """Test organizing an entire project directory."""

    def test_organize_testdata_dir(self, testdata_dir, mock_messages, tmp_path):
        from rc_organize_ras_project import OrganizeRASProject

        output_gdb = str(tmp_path / "TestdataOrganized.gdb")
        tool = OrganizeRASProject()
        params = _build_organize_params(
            testdata_dir, output_gdb,
            include_1d=True, include_2d=True, include_results=True,
        )
        tool.execute(params, mock_messages)

        assert not mock_messages.has_error(), f"Errors: {mock_messages.errors}"


class TestOrganizeLoadToMap:
    """Verify load_to_map parameter controls map loading behavior."""

    def test_load_to_map_false_skips_map(self, hdf_2d_unsteady, mock_messages, tmp_path):
        from rc_organize_ras_project import OrganizeRASProject

        output_gdb = str(tmp_path / "NoMap.gdb")
        tool = OrganizeRASProject()
        params = _build_organize_params(
            hdf_2d_unsteady, output_gdb, load_to_map=False,
        )
        tool.execute(params, mock_messages)

        # Should complete without "Adding results to map..." message
        assert not mock_messages.has_error(), f"Errors: {mock_messages.errors}"
        map_messages = [m for m in mock_messages.messages if "Adding results to map" in m]
        assert len(map_messages) == 0, "Should not attempt map loading when load_to_map=False"


class TestOrganizeNoFiles:
    """Test error handling when no HDF files exist."""

    def test_empty_directory(self, mock_messages, tmp_path):
        from rc_organize_ras_project import OrganizeRASProject

        empty_dir = str(tmp_path / "empty")
        os.makedirs(empty_dir, exist_ok=True)
        output_gdb = str(tmp_path / "Empty.gdb")
        tool = OrganizeRASProject()
        params = _build_organize_params(empty_dir, output_gdb)
        tool.execute(params, mock_messages)

        assert mock_messages.has_error(), "Should report error for empty directory"

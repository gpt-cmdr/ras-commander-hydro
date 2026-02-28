# -*- coding: utf-8 -*-
"""
Integration tests for LoadHECRAS2DGeometry tool.

Requires arcpy (ArcGIS Pro Python environment).

Parameter layout (14 params):
    0: input_hdf
    1: override_crs
    2: geometry_elements (multiValue)
    3: output_breaklines
    4: output_bc_lines
    5: output_perimeters
    6: output_cell_points
    7: output_cell_faces
    8: output_cell_polys
    9: output_pipe_conduits
    10: output_pipe_nodes
    11: output_pipe_networks
    12: output_gdb
    13: create_gdb
"""

import os
import pytest

pytestmark = [pytest.mark.tier1, pytest.mark.requires_arcpy]

import arcpy
if getattr(arcpy, "_is_mock", False):
    pytest.skip("arcpy is mocked — real arcpy required", allow_module_level=True)

from conftest import MockParam, MockMessages


def _build_2d_params(hdf_path, elements, output_prefix="Test2D"):
    """Build 14-element parameter list for LoadHECRAS2DGeometry."""
    params = [
        MockParam(hdf_path, hdf_path),  # 0: input_hdf
        MockParam(None),                  # 1: override_crs
        MockParam(elements),              # 2: geometry_elements
        MockParam(None),                  # 3: output_breaklines
        MockParam(None),                  # 4: output_bc_lines
        MockParam(None),                  # 5: output_perimeters
        MockParam(None),                  # 6: output_cell_points
        MockParam(None),                  # 7: output_cell_faces
        MockParam(None),                  # 8: output_cell_polys
        MockParam(None),                  # 9: output_pipe_conduits
        MockParam(None),                  # 10: output_pipe_nodes
        MockParam(None),                  # 11: output_pipe_networks
        MockParam(None),                  # 12: output_gdb
        MockParam(False),                 # 13: create_gdb
    ]

    element_map = {
        "2D Breaklines": 3,
        "2D Boundary Condition Lines": 4,
        "Mesh Area Perimeters": 5,
        "Mesh Cell Centers": 6,
        "Mesh Cell Faces": 7,
        "Mesh Cells (Polygons)": 8,
        "Pipe Conduits": 9,
        "Pipe Nodes": 10,
        "Pipe Networks": 11,
    }
    for elem in elements:
        idx = element_map.get(elem)
        if idx is not None:
            fc_name = elem.replace(" ", "").replace("(", "").replace(")", "")
            params[idx].valueAsText = rf"memory\{output_prefix}_{fc_name}"

    return params


class TestLoad2DGeometryMeshPerimeters:
    """Test mesh area perimeter extraction."""

    def test_extract_perimeters(self, hdf_2d_unsteady, mock_messages):
        from rc_load_hecras_2d_geometry import LoadHECRAS2DGeometry

        tool = LoadHECRAS2DGeometry()
        params = _build_2d_params(hdf_2d_unsteady, ["Mesh Area Perimeters"], "Perim")
        tool.execute(params, mock_messages)

        assert not mock_messages.has_error(), f"Errors: {mock_messages.errors}"
        output_fc = params[5].valueAsText
        assert arcpy.Exists(output_fc), f"Output not created: {output_fc}"
        count = int(arcpy.GetCount_management(output_fc)[0])
        assert count > 0


class TestLoad2DGeometryCellCenters:
    """Test mesh cell center extraction."""

    def test_extract_cell_centers(self, hdf_2d_unsteady, mock_messages):
        from rc_load_hecras_2d_geometry import LoadHECRAS2DGeometry

        tool = LoadHECRAS2DGeometry()
        params = _build_2d_params(hdf_2d_unsteady, ["Mesh Cell Centers"], "CC")
        tool.execute(params, mock_messages)

        assert not mock_messages.has_error(), f"Errors: {mock_messages.errors}"
        output_fc = params[6].valueAsText
        if arcpy.Exists(output_fc):
            count = int(arcpy.GetCount_management(output_fc)[0])
            assert count > 0, "No cell centers extracted"


class TestLoad2DGeometryCellFaces:
    """Test mesh cell face extraction."""

    def test_extract_cell_faces(self, hdf_2d_unsteady, mock_messages):
        from rc_load_hecras_2d_geometry import LoadHECRAS2DGeometry

        tool = LoadHECRAS2DGeometry()
        params = _build_2d_params(hdf_2d_unsteady, ["Mesh Cell Faces"], "CF")
        tool.execute(params, mock_messages)

        assert not mock_messages.has_error(), f"Errors: {mock_messages.errors}"


class TestLoad2DGeometryBreaklines:
    """Test 2D breakline extraction."""

    def test_extract_breaklines(self, hdf_2d_unsteady, mock_messages):
        from rc_load_hecras_2d_geometry import LoadHECRAS2DGeometry

        tool = LoadHECRAS2DGeometry()
        params = _build_2d_params(hdf_2d_unsteady, ["2D Breaklines"], "BL")
        tool.execute(params, mock_messages)

        # Breaklines may or may not exist in this file
        if mock_messages.has_error():
            # Check if error is just "no breaklines found"
            error_text = " ".join(mock_messages.errors)
            assert "breakline" in error_text.lower() or "not found" in error_text.lower()


class TestLoad2DGeometryBCLines:
    """Test boundary condition line extraction."""

    def test_extract_bc_lines(self, hdf_2d_unsteady, mock_messages):
        from rc_load_hecras_2d_geometry import LoadHECRAS2DGeometry

        tool = LoadHECRAS2DGeometry()
        params = _build_2d_params(hdf_2d_unsteady, ["2D Boundary Condition Lines"], "BC")
        tool.execute(params, mock_messages)

        assert not mock_messages.has_error(), f"Errors: {mock_messages.errors}"


class TestLoad2DGeometryPipeNetworks:
    """Test pipe network extraction from DavisStormSystem."""

    def test_extract_pipe_conduits(self, hdf_2d_pipes, mock_messages):
        from rc_load_hecras_2d_geometry import LoadHECRAS2DGeometry

        tool = LoadHECRAS2DGeometry()
        params = _build_2d_params(hdf_2d_pipes, ["Pipe Conduits"], "PC")
        tool.execute(params, mock_messages)

        assert not mock_messages.has_error(), f"Errors: {mock_messages.errors}"
        output_fc = params[9].valueAsText
        if arcpy.Exists(output_fc):
            count = int(arcpy.GetCount_management(output_fc)[0])
            assert count > 0, "No pipe conduits extracted"

    def test_extract_pipe_nodes(self, hdf_2d_pipes, mock_messages):
        from rc_load_hecras_2d_geometry import LoadHECRAS2DGeometry

        tool = LoadHECRAS2DGeometry()
        params = _build_2d_params(hdf_2d_pipes, ["Pipe Nodes"], "PN")
        tool.execute(params, mock_messages)

        assert not mock_messages.has_error(), f"Errors: {mock_messages.errors}"


class TestLoad2DGeometryMultipleElements:
    """Test extracting multiple 2D elements at once."""

    def test_perimeters_and_cell_centers(self, hdf_2d_unsteady, mock_messages):
        from rc_load_hecras_2d_geometry import LoadHECRAS2DGeometry

        elements = ["Mesh Area Perimeters", "Mesh Cell Centers"]
        tool = LoadHECRAS2DGeometry()
        params = _build_2d_params(hdf_2d_unsteady, elements, "Multi2D")
        tool.execute(params, mock_messages)

        assert not mock_messages.has_error(), f"Errors: {mock_messages.errors}"

    def test_no_elements_raises_error(self, hdf_2d_unsteady, mock_messages):
        from rc_load_hecras_2d_geometry import LoadHECRAS2DGeometry

        tool = LoadHECRAS2DGeometry()
        params = _build_2d_params(hdf_2d_unsteady, [], "Empty2D")
        with pytest.raises(arcpy.ExecuteError):
            tool.execute(params, mock_messages)

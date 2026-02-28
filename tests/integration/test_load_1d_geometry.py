# -*- coding: utf-8 -*-
"""
Integration tests for LoadHECRAS1DGeometry tool.

Requires arcpy (ArcGIS Pro Python environment).
Outputs go to in-memory workspace for speed and cleanup-free testing.

Parameter layout (10 params):
    0: input_hdf (DEFile)
    1: override_crs (GPSpatialReference, optional)
    2: geometry_elements (GPString, multiValue)
    3: output_cross_sections
    4: output_centerlines
    5: output_bank_lines
    6: output_edge_lines
    7: output_structures
    8: output_gdb (optional)
    9: create_gdb (GPBoolean)
"""

import os
import pytest

pytestmark = [pytest.mark.tier1, pytest.mark.requires_arcpy]

import arcpy
if getattr(arcpy, "_is_mock", False):
    pytest.skip("arcpy is mocked — real arcpy required", allow_module_level=True)

from conftest import MockParam, MockMessages


def _build_1d_params(hdf_path, elements, output_prefix="Test"):
    """Build 10-element parameter list for LoadHECRAS1DGeometry."""
    params = [
        MockParam(hdf_path, hdf_path),      # 0: input_hdf
        MockParam(None),                      # 1: override_crs
        MockParam(elements),                  # 2: geometry_elements
        MockParam(None),                      # 3: output_cross_sections
        MockParam(None),                      # 4: output_centerlines
        MockParam(None),                      # 5: output_bank_lines
        MockParam(None),                      # 6: output_edge_lines
        MockParam(None),                      # 7: output_structures
        MockParam(None),                      # 8: output_gdb
        MockParam(False),                     # 9: create_gdb
    ]

    # Set output paths to in-memory workspace for selected elements
    element_map = {
        "Cross Sections": 3,
        "River Centerlines": 4,
        "Bank Lines": 5,
        "Edge Lines": 6,
        "1D Structures": 7,
    }
    for elem in elements:
        idx = element_map.get(elem)
        if idx is not None:
            fc_name = elem.replace(" ", "")
            params[idx].valueAsText = rf"memory\{output_prefix}_{fc_name}"

    return params


class TestLoad1DGeometryCrossSections:
    """Test cross section extraction from 1D HDF files."""

    def test_extract_cross_sections(self, hdf_1d_unsteady, mock_messages):
        from rc_load_hecras_1d_geometry import LoadHECRAS1DGeometry

        tool = LoadHECRAS1DGeometry()
        params = _build_1d_params(hdf_1d_unsteady, ["Cross Sections"])
        tool.execute(params, mock_messages)

        assert not mock_messages.has_error(), f"Errors: {mock_messages.errors}"
        # Verify output was created
        output_fc = params[3].valueAsText
        assert arcpy.Exists(output_fc), f"Output FC not created: {output_fc}"
        count = int(arcpy.GetCount_management(output_fc)[0])
        assert count > 0, "No cross sections extracted"

    def test_cross_sections_have_fields(self, hdf_1d_unsteady, mock_messages):
        from rc_load_hecras_1d_geometry import LoadHECRAS1DGeometry

        tool = LoadHECRAS1DGeometry()
        params = _build_1d_params(hdf_1d_unsteady, ["Cross Sections"], "Fields")
        tool.execute(params, mock_messages)

        output_fc = params[3].valueAsText
        if arcpy.Exists(output_fc):
            field_names = [f.name for f in arcpy.ListFields(output_fc)]
            # Should have at minimum River and Reach fields
            assert "River" in field_names or "xs_id" in field_names

    def test_steady_cross_sections(self, hdf_1d_steady, mock_messages):
        from rc_load_hecras_1d_geometry import LoadHECRAS1DGeometry

        tool = LoadHECRAS1DGeometry()
        params = _build_1d_params(hdf_1d_steady, ["Cross Sections"], "Steady")
        tool.execute(params, mock_messages)

        assert not mock_messages.has_error(), f"Errors: {mock_messages.errors}"


class TestLoad1DGeometryCenterlines:
    """Test river centerline extraction."""

    def test_extract_centerlines(self, hdf_1d_unsteady, mock_messages):
        from rc_load_hecras_1d_geometry import LoadHECRAS1DGeometry

        tool = LoadHECRAS1DGeometry()
        params = _build_1d_params(hdf_1d_unsteady, ["River Centerlines"], "CL")
        tool.execute(params, mock_messages)

        assert not mock_messages.has_error(), f"Errors: {mock_messages.errors}"
        output_fc = params[4].valueAsText
        if arcpy.Exists(output_fc):
            count = int(arcpy.GetCount_management(output_fc)[0])
            assert count > 0


class TestLoad1DGeometryMultipleElements:
    """Test extracting multiple elements at once."""

    def test_all_elements(self, hdf_1d_unsteady, mock_messages):
        from rc_load_hecras_1d_geometry import LoadHECRAS1DGeometry

        elements = ["Cross Sections", "River Centerlines", "Bank Lines", "Edge Lines"]
        tool = LoadHECRAS1DGeometry()
        params = _build_1d_params(hdf_1d_unsteady, elements, "Multi")
        tool.execute(params, mock_messages)

        assert not mock_messages.has_error(), f"Errors: {mock_messages.errors}"

    def test_no_elements_raises_error(self, hdf_1d_unsteady, mock_messages):
        from rc_load_hecras_1d_geometry import LoadHECRAS1DGeometry

        tool = LoadHECRAS1DGeometry()
        params = _build_1d_params(hdf_1d_unsteady, [], "Empty")
        with pytest.raises(arcpy.ExecuteError):
            tool.execute(params, mock_messages)

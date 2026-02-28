# -*- coding: utf-8 -*-
"""
Integration tests for LoadHECRAS2DResults tool.

Requires arcpy (ArcGIS Pro Python environment).

Parameter layout (7 params):
    0: input_hdf
    1: override_crs
    2: results_elements (multiValue)
    3: output_max_wse
    4: output_max_face_vel
    5: output_gdb
    6: create_gdb
"""

import os
import pytest

pytestmark = [pytest.mark.tier1, pytest.mark.requires_arcpy]

import arcpy
if getattr(arcpy, "_is_mock", False):
    pytest.skip("arcpy is mocked — real arcpy required", allow_module_level=True)

from conftest import MockParam, MockMessages


def _build_results_params(hdf_path, elements, output_prefix="TestRes"):
    """Build 7-element parameter list for LoadHECRAS2DResults."""
    params = [
        MockParam(hdf_path, hdf_path),  # 0: input_hdf
        MockParam(None),                  # 1: override_crs
        MockParam(elements),              # 2: results_elements
        MockParam(None),                  # 3: output_max_wse
        MockParam(None),                  # 4: output_max_face_vel
        MockParam(None),                  # 5: output_gdb
        MockParam(False),                 # 6: create_gdb
    ]

    element_map = {
        "Max WSE at Cell Centers": 3,
        "Max Vel at Cell Faces": 4,
    }
    for elem in elements:
        idx = element_map.get(elem)
        if idx is not None:
            fc_name = elem.replace(" ", "")
            params[idx].valueAsText = rf"memory\{output_prefix}_{fc_name}"

    return params


class TestLoadMaxWSE:
    """Test maximum water surface elevation extraction."""

    def test_extract_max_wse(self, hdf_2d_unsteady, mock_messages):
        from rc_load_hecras_2d_results import LoadHECRAS2DResults

        tool = LoadHECRAS2DResults()
        params = _build_results_params(
            hdf_2d_unsteady, ["Max WSE at Cell Centers"], "WSE"
        )
        tool.execute(params, mock_messages)

        assert not mock_messages.has_error(), f"Errors: {mock_messages.errors}"
        output_fc = params[3].valueAsText
        assert arcpy.Exists(output_fc), f"Output not created: {output_fc}"
        count = int(arcpy.GetCount_management(output_fc)[0])
        assert count > 0, "No WSE points extracted"

    def test_wse_has_expected_fields(self, hdf_2d_unsteady, mock_messages):
        from rc_load_hecras_2d_results import LoadHECRAS2DResults

        tool = LoadHECRAS2DResults()
        params = _build_results_params(
            hdf_2d_unsteady, ["Max WSE at Cell Centers"], "WSEFields"
        )
        tool.execute(params, mock_messages)

        output_fc = params[3].valueAsText
        if arcpy.Exists(output_fc):
            field_names = [f.name for f in arcpy.ListFields(output_fc)]
            # Should have WSE-related fields
            assert len(field_names) > 2  # More than just OID + Shape


class TestLoadMaxVelocity:
    """Test maximum velocity extraction."""

    def test_extract_max_velocity(self, hdf_2d_unsteady, mock_messages):
        from rc_load_hecras_2d_results import LoadHECRAS2DResults

        tool = LoadHECRAS2DResults()
        params = _build_results_params(
            hdf_2d_unsteady, ["Max Vel at Cell Faces"], "Vel"
        )
        tool.execute(params, mock_messages)

        assert not mock_messages.has_error(), f"Errors: {mock_messages.errors}"


class TestLoadBothResults:
    """Test extracting both WSE and velocity."""

    def test_extract_both(self, hdf_2d_unsteady, mock_messages):
        from rc_load_hecras_2d_results import LoadHECRAS2DResults

        elements = ["Max WSE at Cell Centers", "Max Vel at Cell Faces"]
        tool = LoadHECRAS2DResults()
        params = _build_results_params(hdf_2d_unsteady, elements, "Both")
        tool.execute(params, mock_messages)

        assert not mock_messages.has_error(), f"Errors: {mock_messages.errors}"

    def test_no_elements_raises_error(self, hdf_2d_unsteady, mock_messages):
        from rc_load_hecras_2d_results import LoadHECRAS2DResults

        tool = LoadHECRAS2DResults()
        params = _build_results_params(hdf_2d_unsteady, [], "EmptyRes")
        with pytest.raises(arcpy.ExecuteError):
            tool.execute(params, mock_messages)


class TestResultsOn1DFile:
    """Test results tool on 1D file (which has no 2D results)."""

    def test_1d_file_produces_no_2d_results(self, hdf_1d_unsteady, mock_messages):
        """1D file has Plan Info (so has_results=True) but no 2D flow area results.

        The tool proceeds past the has_results check but finds no 2D data to
        extract. It either raises ExecuteError or produces warnings/no output.
        """
        from rc_load_hecras_2d_results import LoadHECRAS2DResults

        tool = LoadHECRAS2DResults()
        params = _build_results_params(
            hdf_1d_unsteady, ["Max WSE at Cell Centers"], "NoRes"
        )
        try:
            tool.execute(params, mock_messages)
        except arcpy.ExecuteError:
            pass  # Expected — tool may raise on missing 2D data

        # Regardless of whether it raised, it should not have created the output
        output_fc = params[3].valueAsText
        if output_fc:
            assert not arcpy.Exists(output_fc) or mock_messages.has_error() or mock_messages.has_warning()

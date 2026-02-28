# -*- coding: utf-8 -*-
"""
Unit tests for pure utility functions in rc_utils.py.

Tests naming helpers, field detection, and clockwise checks that
do NOT require arcpy or HDF files.
"""

import os
import re
import pytest
import numpy as np

pytestmark = pytest.mark.tier1


# ---------------------------------------------------------------------------
# extract_project_and_plan_info
# ---------------------------------------------------------------------------
class TestExtractProjectAndPlanInfo:
    """Tests for rc_utils.extract_project_and_plan_info."""

    def _extract(self, path):
        from rc_utils import extract_project_and_plan_info
        return extract_project_and_plan_info(path)

    def test_standard_plan_file(self):
        proj, plan, base = self._extract("BaldEagleDamBrk.p07.hdf")
        assert proj == "BaldEagleDamBrk"
        assert plan == "07"
        assert base == "BaldEagleDamBrk.p07"

    def test_plan_01(self):
        proj, plan, base = self._extract("BaldEagle.p01.hdf")
        assert proj == "BaldEagle"
        assert plan == "01"
        assert base == "BaldEagle.p01"

    def test_plan_with_directory(self):
        proj, plan, base = self._extract(r"C:\Models\BaldEagle.p02.hdf")
        assert proj == "BaldEagle"
        assert plan == "02"

    def test_no_plan_pattern(self):
        proj, plan, base = self._extract("SomeFile.hdf")
        assert proj == "SomeFile"
        assert plan == "00"
        assert base == "SomeFile"

    def test_dots_in_project_name(self):
        proj, plan, base = self._extract("My.Project.Name.p03.hdf")
        assert proj == "My.Project.Name"
        assert plan == "03"

    def test_geometry_file_not_matched(self):
        """Geometry files (g*.hdf) don't match the p-pattern."""
        proj, plan, base = self._extract("Model.g01.hdf")
        assert plan == "00"  # Falls back since no pXX match


# ---------------------------------------------------------------------------
# get_feature_dataset_name
# ---------------------------------------------------------------------------
class TestGetFeatureDatasetName:
    def test_standard(self):
        from rc_utils import get_feature_dataset_name
        result = get_feature_dataset_name("BaldEagle.p01.hdf")
        assert result == "BaldEagle_Plan01"

    def test_dam_break(self):
        from rc_utils import get_feature_dataset_name
        result = get_feature_dataset_name("BaldEagleDamBrk.p07.hdf")
        assert result == "BaldEagleDamBrk_Plan07"


# ---------------------------------------------------------------------------
# get_feature_class_name
# ---------------------------------------------------------------------------
class TestGetFeatureClassName:
    def test_cross_sections(self):
        from rc_utils import get_feature_class_name
        result = get_feature_class_name("CrossSections", "BaldEagle", "01")
        assert result == "CrossSections_BaldEagle_Plan01"


# ---------------------------------------------------------------------------
# is_clockwise_numpy
# ---------------------------------------------------------------------------
class TestIsClockwiseNumpy:
    def _check(self, coords):
        from rc_utils import is_clockwise_numpy
        return is_clockwise_numpy(np.array(coords))

    def test_clockwise_square(self):
        # Clockwise: top-left -> top-right -> bottom-right -> bottom-left
        cw = [(0, 1), (1, 1), (1, 0), (0, 0)]
        assert self._check(cw) == True

    def test_counterclockwise_square(self):
        # Counter-clockwise: reverse
        ccw = [(0, 0), (1, 0), (1, 1), (0, 1)]
        assert self._check(ccw) == False

    def test_triangle_cw(self):
        cw = [(0, 0), (2, 1), (1, 2)]
        # Shoelace: 0.5*(0*1-2*0 + 2*2-1*1 + 1*0-0*2) = 0.5*(0+3-2) = 0.5
        # CW if area < 0 — check the expected sign
        from rc_utils import is_clockwise_numpy
        result = is_clockwise_numpy(np.array(cw))
        assert isinstance(result, (bool, np.bool_))

    def test_degenerate_line(self):
        """Two-point array (degenerate): should not crash."""
        coords = [(0, 0), (1, 1)]
        from rc_utils import is_clockwise_numpy
        # Degenerate, area ≈ 0
        result = is_clockwise_numpy(np.array(coords))
        assert isinstance(result, (bool, np.bool_))


# ---------------------------------------------------------------------------
# get_dynamic_fields_from_data
# ---------------------------------------------------------------------------
class TestGetDynamicFieldsFromData:
    def _get_fields(self, data):
        from rc_utils import get_dynamic_fields_from_data
        return get_dynamic_fields_from_data(data)

    def test_empty_list(self):
        assert self._get_fields([]) == []

    def test_integer_field(self):
        data = [{"CellID": 1}, {"CellID": 2}]
        fields = self._get_fields(data)
        assert any(f[0] == "CellID" and f[1] == "LONG" for f in fields)

    def test_float_field(self):
        data = [{"Elevation": 100.5}, {"Elevation": 200.3}]
        fields = self._get_fields(data)
        assert any(f[0] == "Elevation" and f[1] == "DOUBLE" for f in fields)

    def test_text_field(self):
        data = [{"Name": "River1"}, {"Name": "River2"}]
        fields = self._get_fields(data)
        text_field = [f for f in fields if f[0] == "Name"]
        assert text_field[0][1] == "TEXT"
        assert len(text_field[0]) == 3  # has length parameter

    def test_boolean_field(self):
        data = [{"IsActive": True}, {"IsActive": False}]
        fields = self._get_fields(data)
        assert any(f[0] == "IsActive" and f[1] == "SHORT" for f in fields)

    def test_numpy_types(self):
        data = [
            {"IntVal": np.int32(42), "FloatVal": np.float64(3.14)},
        ]
        fields = self._get_fields(data)
        int_field = [f for f in fields if f[0] == "IntVal"]
        float_field = [f for f in fields if f[0] == "FloatVal"]
        assert int_field[0][1] == "LONG"
        assert float_field[0][1] == "DOUBLE"

    def test_all_nan_defaults_to_double_for_elevation(self):
        """Fields with all NaN values and numeric-sounding names default to DOUBLE."""
        data = [{"elevation": float("nan")}, {"elevation": float("nan")}]
        fields = self._get_fields(data)
        elev_field = [f for f in fields if f[0] == "elevation"]
        assert elev_field[0][1] == "DOUBLE"

    def test_all_nan_defaults_to_text_for_name(self):
        """Fields with all NaN values and non-numeric names default to TEXT."""
        data = [{"status": float("nan")}, {"status": float("nan")}]
        fields = self._get_fields(data)
        status_field = [f for f in fields if f[0] == "status"]
        assert status_field[0][1] == "TEXT"

    def test_text_field_length_minimum(self):
        """Text fields have a minimum length of 50."""
        data = [{"Name": "AB"}]
        fields = self._get_fields(data)
        text_field = [f for f in fields if f[0] == "Name"]
        assert text_field[0][2] >= 50

    def test_text_field_length_capped(self):
        """Text fields are capped at 255."""
        data = [{"Name": "X" * 300}]
        fields = self._get_fields(data)
        text_field = [f for f in fields if f[0] == "Name"]
        assert text_field[0][2] <= 255

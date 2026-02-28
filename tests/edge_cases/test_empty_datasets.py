# -*- coding: utf-8 -*-
"""
Edge case tests for empty or degenerate HDF datasets.

Covers:
    - Zero-length arrays in HDF datasets
    - All-NaN fields in attribute tables
    - Missing expected datasets within existing groups
"""

import os
import pytest
import h5py
import numpy as np

pytestmark = [pytest.mark.tier1]


class TestZeroLengthDatasets:
    """Test handling of zero-length HDF datasets."""

    @pytest.fixture
    def hdf_with_empty_datasets(self, tmp_path):
        """Create HDF file with empty datasets inside expected group structure."""
        path = str(tmp_path / "empty_datasets.hdf")
        with h5py.File(path, "w") as f:
            # Create geometry group with empty cross sections
            xs_group = f.create_group("Geometry/Cross Sections")
            # Empty polyline points (0 rows, 2 columns)
            xs_group.create_dataset("Polyline Points", shape=(0, 2), dtype="float64")
            # Empty attributes structured array
            dt = np.dtype([("Name", "S64"), ("River", "S64"), ("Reach", "S64")])
            xs_group.create_dataset("Attributes", shape=(0,), dtype=dt)

            # Create 2D flow area group with empty cells
            fa_group = f.create_group("Geometry/2D Flow Areas")
            attr_dt = np.dtype([("Name", "S64")])
            attrs = np.array([(b"TestMesh",)], dtype=attr_dt)
            fa_group.create_dataset("Attributes", data=attrs)

            mesh_group = f.create_group("Geometry/2D Flow Areas/TestMesh")
            mesh_group.create_dataset("Cells Center Coordinate", shape=(0, 2), dtype="float64")
            mesh_group.create_dataset("Faces FacePoint Indexes", shape=(0, 2), dtype="int32")
        return path

    def test_empty_cross_section_attributes(self, hdf_with_empty_datasets):
        with h5py.File(hdf_with_empty_datasets, "r") as f:
            attrs = f["Geometry/Cross Sections/Attributes"][()]
            assert len(attrs) == 0

    def test_empty_polyline_points(self, hdf_with_empty_datasets):
        with h5py.File(hdf_with_empty_datasets, "r") as f:
            pts = f["Geometry/Cross Sections/Polyline Points"][()]
            assert pts.shape[0] == 0

    def test_empty_cell_centers(self, hdf_with_empty_datasets):
        with h5py.File(hdf_with_empty_datasets, "r") as f:
            centers = f["Geometry/2D Flow Areas/TestMesh/Cells Center Coordinate"][()]
            assert len(centers) == 0

    def test_check_contents_with_empty_datasets(self, hdf_with_empty_datasets):
        from rc_organize_ras_project import OrganizeRASProject

        tool = OrganizeRASProject()
        with h5py.File(hdf_with_empty_datasets, "r") as f:
            contents = tool._check_hdf_contents(f)
        # Cross sections group exists with Attributes, so has_1d should be True
        assert contents["has_1d"] is True
        # 2D flow areas exist
        assert contents["has_2d"] is True


class TestAllNaNFields:
    """Test get_dynamic_fields_from_data with all-NaN values."""

    def test_all_nan_records(self):
        from rc_utils import get_dynamic_fields_from_data

        data = [
            {"elevation": float("nan"), "area": float("nan")},
            {"elevation": float("nan"), "area": float("nan")},
        ]
        fields = get_dynamic_fields_from_data(data)
        assert len(fields) > 0
        # Both should default to DOUBLE since they have numeric keywords
        for name, ftype, *_ in fields:
            assert ftype == "DOUBLE"

    def test_mixed_nan_and_values(self):
        from rc_utils import get_dynamic_fields_from_data

        data = [
            {"elevation": float("nan"), "name": "River1"},
            {"elevation": 100.5, "name": float("nan")},
        ]
        fields = get_dynamic_fields_from_data(data)
        elev_field = [f for f in fields if f[0] == "elevation"]
        name_field = [f for f in fields if f[0] == "name"]
        assert elev_field[0][1] == "DOUBLE"
        assert name_field[0][1] == "TEXT"

    def test_numpy_nan_values(self):
        from rc_utils import get_dynamic_fields_from_data

        data = [
            {"value": np.float64("nan")},
            {"value": np.float64(42.0)},
        ]
        fields = get_dynamic_fields_from_data(data)
        val_field = [f for f in fields if f[0] == "value"]
        assert val_field[0][1] == "DOUBLE"


class TestMissingExpectedDatasets:
    """Test HDF files where parent groups exist but expected child datasets don't."""

    @pytest.fixture
    def hdf_missing_polylines(self, tmp_path):
        """Cross Sections group exists but Polyline Points dataset is missing."""
        path = str(tmp_path / "missing_polylines.hdf")
        with h5py.File(path, "w") as f:
            xs_group = f.create_group("Geometry/Cross Sections")
            dt = np.dtype([("Name", "S64"), ("River", "S64"), ("Reach", "S64")])
            attrs = np.array([(b"XS1", b"River1", b"Reach1")], dtype=dt)
            xs_group.create_dataset("Attributes", data=attrs)
            # Intentionally NOT creating Polyline Points
        return path

    def test_attributes_without_polylines(self, hdf_missing_polylines):
        with h5py.File(hdf_missing_polylines, "r") as f:
            assert "Geometry/Cross Sections/Attributes" in f
            assert "Geometry/Cross Sections/Polyline Points" not in f

    @pytest.fixture
    def hdf_missing_face_indexes(self, tmp_path):
        """2D Flow Area with cell centers but no face indexes."""
        path = str(tmp_path / "missing_faces.hdf")
        with h5py.File(path, "w") as f:
            fa_group = f.create_group("Geometry/2D Flow Areas")
            attr_dt = np.dtype([("Name", "S64")])
            attrs = np.array([(b"TestMesh",)], dtype=attr_dt)
            fa_group.create_dataset("Attributes", data=attrs)

            mesh_group = f.create_group("Geometry/2D Flow Areas/TestMesh")
            centers = np.array([[1.0, 2.0], [3.0, 4.0]])
            mesh_group.create_dataset("Cells Center Coordinate", data=centers)
            # Intentionally NOT creating Faces FacePoint Indexes
        return path

    def test_cell_centers_without_faces(self, hdf_missing_face_indexes):
        with h5py.File(hdf_missing_face_indexes, "r") as f:
            assert "Geometry/2D Flow Areas/TestMesh/Cells Center Coordinate" in f
            assert "Geometry/2D Flow Areas/TestMesh/Faces FacePoint Indexes" not in f


class TestCacheHDFMetadataEdgeCases:
    """Test cache_hdf_metadata with unusual HDF structures."""

    def test_cache_no_results(self, hdf_1d_steady):
        """Steady file may not have simulation start time."""
        from rc_utils import cache_hdf_metadata

        with h5py.File(hdf_1d_steady, "r") as f:
            cache = cache_hdf_metadata(f)
        # Should not crash — results may or may not be present
        assert isinstance(cache, dict)
        assert "has_results" in cache

    def test_cache_with_2d_results(self, hdf_2d_unsteady):
        from rc_utils import cache_hdf_metadata

        with h5py.File(hdf_2d_unsteady, "r") as f:
            cache = cache_hdf_metadata(f)
        assert cache["has_results"] is True
        assert cache["simulation_start_time"] is not None
        assert len(cache["mesh_names"]) > 0

    def test_cache_with_pipes(self, hdf_2d_pipes):
        from rc_utils import cache_hdf_metadata

        with h5py.File(hdf_2d_pipes, "r") as f:
            cache = cache_hdf_metadata(f)
        assert cache["has_pipe_conduits"] is True
        assert cache["has_pipe_nodes"] is True

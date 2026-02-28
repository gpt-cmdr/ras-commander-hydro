# -*- coding: utf-8 -*-
"""
Unit tests for HDF structure validation using h5py directly.

These tests verify HDF file structure and metadata without creating
arcpy geometry objects, so they can run without the arcpy license
in scenarios where only h5py is available.
"""

import os
import pytest
import h5py
import numpy as np

pytestmark = [pytest.mark.tier1]


# ---------------------------------------------------------------------------
# HDF Structure — 1D Geometry
# ---------------------------------------------------------------------------
class TestHDF1DStructure:
    """Validate expected HDF groups/datasets in 1D plan files."""

    def test_has_geometry_group(self, hdf_1d_unsteady):
        with h5py.File(hdf_1d_unsteady, "r") as f:
            assert "Geometry" in f

    def test_has_cross_sections(self, hdf_1d_unsteady):
        with h5py.File(hdf_1d_unsteady, "r") as f:
            assert "Geometry/Cross Sections" in f

    def test_cross_section_attributes_exist(self, hdf_1d_unsteady):
        with h5py.File(hdf_1d_unsteady, "r") as f:
            assert "Geometry/Cross Sections/Attributes" in f

    def test_has_river_centerlines(self, hdf_1d_unsteady):
        with h5py.File(hdf_1d_unsteady, "r") as f:
            assert "Geometry/River Centerlines" in f

    def test_cross_section_polylines(self, hdf_1d_unsteady):
        with h5py.File(hdf_1d_unsteady, "r") as f:
            path = "Geometry/Cross Sections/Polyline Points"
            if path in f:
                pts = f[path][()]
                assert pts.shape[1] >= 2  # at least x, y

    def test_steady_also_has_geometry(self, hdf_1d_steady):
        with h5py.File(hdf_1d_steady, "r") as f:
            assert "Geometry" in f
            assert "Geometry/Cross Sections" in f

    def test_plan_data_exists_unsteady(self, hdf_1d_unsteady):
        with h5py.File(hdf_1d_unsteady, "r") as f:
            assert "Plan Data" in f

    def test_plan_info_exists_unsteady(self, hdf_1d_unsteady):
        with h5py.File(hdf_1d_unsteady, "r") as f:
            assert "Plan Data/Plan Information" in f


# ---------------------------------------------------------------------------
# HDF Structure — 2D Geometry
# ---------------------------------------------------------------------------
class TestHDF2DStructure:
    """Validate expected HDF groups/datasets in 2D plan files."""

    def test_has_2d_flow_areas(self, hdf_2d_unsteady):
        with h5py.File(hdf_2d_unsteady, "r") as f:
            assert "Geometry/2D Flow Areas" in f

    def test_2d_flow_areas_have_attributes(self, hdf_2d_unsteady):
        with h5py.File(hdf_2d_unsteady, "r") as f:
            assert "Geometry/2D Flow Areas/Attributes" in f

    def test_mesh_has_cell_centers(self, hdf_2d_unsteady):
        with h5py.File(hdf_2d_unsteady, "r") as f:
            fa_path = "Geometry/2D Flow Areas"
            attrs = f[f"{fa_path}/Attributes"][()]
            mesh_name = attrs["Name"][0].decode("utf-8", "ignore").strip()
            center_path = f"{fa_path}/{mesh_name}/Cells Center Coordinate"
            assert center_path in f
            centers = f[center_path][()]
            assert len(centers) > 0

    def test_mesh_has_face_indexes(self, hdf_2d_unsteady):
        with h5py.File(hdf_2d_unsteady, "r") as f:
            fa_path = "Geometry/2D Flow Areas"
            attrs = f[f"{fa_path}/Attributes"][()]
            mesh_name = attrs["Name"][0].decode("utf-8", "ignore").strip()
            face_path = f"{fa_path}/{mesh_name}/Faces FacePoint Indexes"
            assert face_path in f

    def test_has_results_unsteady(self, hdf_2d_unsteady):
        with h5py.File(hdf_2d_unsteady, "r") as f:
            results_path = "Results/Unsteady/Output/Output Blocks/Base Output/Summary Output/2D Flow Areas"
            assert results_path in f


# ---------------------------------------------------------------------------
# HDF Structure — Pipe Networks
# ---------------------------------------------------------------------------
class TestHDFPipeStructure:
    """Validate pipe network data in DavisStormSystem."""

    def test_has_pipe_conduits(self, hdf_2d_pipes):
        with h5py.File(hdf_2d_pipes, "r") as f:
            assert "Geometry/Pipe Conduits" in f

    def test_has_pipe_nodes(self, hdf_2d_pipes):
        with h5py.File(hdf_2d_pipes, "r") as f:
            assert "Geometry/Pipe Nodes" in f

    def test_pipe_conduits_have_attributes(self, hdf_2d_pipes):
        with h5py.File(hdf_2d_pipes, "r") as f:
            path = "Geometry/Pipe Conduits/Attributes"
            if path in f:
                attrs = f[path][()]
                assert len(attrs) > 0

    def test_pipe_nodes_have_attributes(self, hdf_2d_pipes):
        with h5py.File(hdf_2d_pipes, "r") as f:
            path = "Geometry/Pipe Nodes/Attributes"
            if path in f:
                attrs = f[path][()]
                assert len(attrs) > 0


# ---------------------------------------------------------------------------
# HDF Projection Attribute
# ---------------------------------------------------------------------------
class TestHDFProjection:
    """Verify projection attribute reading from HDF files."""

    def test_projection_attribute_readable(self, hdf_2d_unsteady):
        with h5py.File(hdf_2d_unsteady, "r") as f:
            proj = f.attrs.get("Projection")
            if proj is not None:
                if isinstance(proj, (bytes, np.bytes_)):
                    wkt = proj.decode("utf-8", "ignore")
                else:
                    wkt = str(proj)
                assert len(wkt) > 0

    def test_projection_on_1d_file(self, hdf_1d_unsteady):
        """1D files may or may not have a projection attribute."""
        with h5py.File(hdf_1d_unsteady, "r") as f:
            proj = f.attrs.get("Projection")
            # Just verify no crash - projection may be absent for 1D


# ---------------------------------------------------------------------------
# Simulation Start Time
# ---------------------------------------------------------------------------
class TestSimulationStartTime:
    """Verify simulation start time parsing from HDF metadata."""

    def test_start_time_readable(self, hdf_2d_unsteady):
        with h5py.File(hdf_2d_unsteady, "r") as f:
            plan_info = f.get("Plan Data/Plan Information")
            if plan_info and "Simulation Start Time" in plan_info.attrs:
                time_str = plan_info.attrs["Simulation Start Time"]
                if isinstance(time_str, (bytes, np.bytes_)):
                    decoded = time_str.decode("utf-8", "ignore")
                else:
                    decoded = str(time_str)
                assert len(decoded) > 0

    def test_start_time_parseable(self, hdf_2d_unsteady):
        from datetime import datetime

        with h5py.File(hdf_2d_unsteady, "r") as f:
            plan_info = f.get("Plan Data/Plan Information")
            if plan_info and "Simulation Start Time" in plan_info.attrs:
                time_str = plan_info.attrs["Simulation Start Time"]
                if isinstance(time_str, (bytes, np.bytes_)):
                    decoded = time_str.decode("utf-8", "ignore")
                else:
                    decoded = str(time_str)
                # Should parse as ddMONyyyy HH:MM:SS
                dt = datetime.strptime(decoded, "%d%b%Y %H:%M:%S")
                assert dt.year > 1900

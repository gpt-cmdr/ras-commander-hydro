# -*- coding: utf-8 -*-
"""
LoadHECRAS6xHDFData.py

Tool for loading HEC-RAS 6.x HDF data including 2D geometry elements and results.
This tool's logic is derived from the HdfMesh and HdfBndry classes
in the ras-commander library.
"""

import arcpy
import os
import h5py
import numpy as np
from datetime import datetime, timedelta
from collections import defaultdict

# Import helper functions from utils
from .utils import (
    get_ras_projection_wkt,
    polygonize_arcpy_optimized,
    get_polyline_centroid_vectorized
)


class LoadHECRAS6xHDFData(object):
    """
    Loads one or more 2D geometry element types from a HEC-RAS HDF file.
    This tool's logic is derived from the HdfMesh and HdfBndry classes
    in the ras-commander library.
    """
    def __init__(self):
        self.label = "Load HEC-RAS 6.x HDF Data"
        self.description = "Extracts selected 2D geometry elements and results from a HEC-RAS HDF file into separate feature classes."
        self.canRunInBackground = False
        
        # Geometry elements
        self.BREAKLINES = "2D Breaklines"
        self.BC_LINES = "2D Boundary Condition Lines"
        self.PERIMETERS = "Mesh Area Perimeters"
        self.CELL_POINTS = "Mesh Cell Centers"
        self.CELL_FACES = "Mesh Cell Faces"
        self.CELL_POLYS = "Mesh Cells (Polygons)"
        
        # Pipe network elements
        self.PIPE_CONDUITS = "Pipe Conduits"
        self.PIPE_NODES = "Pipe Nodes"
        self.PIPE_NETWORKS = "Pipe Networks"
        
        # Results elements
        self.MAX_WSE_POINTS = "Max WSE at Cell Centers"
        self.MAX_FACE_VEL_POINTS = "Max Vel at Cell Faces"
        
        # Cache for HDF metadata
        self._hdf_cache = {}

    def getParameterInfo(self):
        # Separate lists for geometry and results
        geometry_elements = [self.BREAKLINES, self.BC_LINES, self.PERIMETERS, self.CELL_POINTS, 
                           self.CELL_FACES, self.CELL_POLYS, self.PIPE_CONDUITS, self.PIPE_NODES, self.PIPE_NETWORKS]
        results_elements = [self.MAX_WSE_POINTS, self.MAX_FACE_VEL_POINTS]

        params = [
            arcpy.Parameter(displayName="Geometry or Plan HDF File", name="input_hdf", datatype="DEFile", 
                          parameterType="Required", direction="Input"),
            arcpy.Parameter(displayName="Override CRS (Optional)", name="override_crs", datatype="GPSpatialReference", 
                          parameterType="Optional", direction="Input"),
            
            # Geometry elements to load
            arcpy.Parameter(displayName="Geometry Elements to Load", name="geometry_elements", datatype="GPString", 
                          parameterType="Optional", direction="Input", multiValue=True, 
                          category="Load 2D Geometry Elements"),
            
            # Results elements to load
            arcpy.Parameter(displayName="Results to Load", name="results_elements", datatype="GPString", 
                          parameterType="Optional", direction="Input", multiValue=True, 
                          category="Load 2D Results"),
            
            # Output parameters
            arcpy.Parameter(displayName="Output 2D Breaklines", name="output_breaklines", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output 2D Boundary Condition Lines", name="output_bc_lines", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output Mesh Area Perimeters", name="output_perimeters", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output Mesh Cell Centers", name="output_cell_points", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output Mesh Cell Faces", name="output_cell_faces", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output Mesh Cells (Polygons)", name="output_cell_polys", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output Max WSE at Cell Centers", name="output_max_wse", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output Max Vel at Cell Faces", name="output_max_face_vel", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output Pipe Conduits", name="output_pipe_conduits", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output Pipe Nodes", name="output_pipe_nodes", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output Pipe Networks", name="output_pipe_networks", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs")
        ]
        
        # Configure HDF file filter
        params[0].filter.list = ["hdf", "g*.hdf", "p*.hdf"]
        
        # Set filters for multi-value parameters
        params[2].filter.type = "ValueList"
        params[2].filter.list = geometry_elements
        params[2].value = [self.PERIMETERS]  # Default selection
        
        params[3].filter.type = "ValueList"
        params[3].filter.list = results_elements
        
        # Set default output paths
        params[4].value = r"memory\Breaklines"
        params[5].value = r"memory\BoundaryConditionLines"
        params[6].value = r"memory\MeshPerimeters"
        params[7].value = r"memory\MeshCellCenters"
        params[8].value = r"memory\MeshCellFaces"
        params[9].value = r"memory\MeshCellPolygons"
        params[10].value = r"memory\MaximumWSE"
        params[11].value = r"memory\MaximumFaceVelocity"
        params[12].value = r"memory\PipeConduits"
        params[13].value = r"memory\PipeNodes"
        params[14].value = r"memory\PipeNetworks"
        
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def _cache_hdf_metadata(self, hdf_file):
        """Pre-cache HDF metadata for faster access."""
        self._hdf_cache = {
            'mesh_names': [],
            'mesh_metadata': {},
            'has_results': False,
            'simulation_start_time': None
        }
        
        # Get mesh names
        flow_areas_path = "Geometry/2D Flow Areas"
        if flow_areas_path in hdf_file and "Attributes" in hdf_file[flow_areas_path]:
            attributes = hdf_file[f"{flow_areas_path}/Attributes"][()]
            self._hdf_cache['mesh_names'] = [n.decode('utf-8', 'ignore').strip() 
                                           for n in attributes["Name"]]
            
            # Cache mesh metadata
            for mesh_name in self._hdf_cache['mesh_names']:
                base_path = f"{flow_areas_path}/{mesh_name}"
                metadata = {}
                
                # Cache dataset sizes
                if f"{base_path}/Cells Center Coordinate" in hdf_file:
                    metadata['cell_count'] = len(hdf_file[f"{base_path}/Cells Center Coordinate"])
                
                if f"{base_path}/Faces FacePoint Indexes" in hdf_file:
                    metadata['face_count'] = len(hdf_file[f"{base_path}/Faces FacePoint Indexes"])
                
                self._hdf_cache['mesh_metadata'][mesh_name] = metadata
        
        # Check for boundary condition lines
        bc_lines_path = "Geometry/Boundary Condition Lines"
        if bc_lines_path in hdf_file:
            self._hdf_cache['has_bc_lines'] = True
            self._hdf_cache['bc_lines_count'] = len(hdf_file[f"{bc_lines_path}/Attributes"][()])
        else:
            self._hdf_cache['has_bc_lines'] = False
        
        # Check for pipe network elements
        pipe_conduits_path = "Geometry/Pipe Conduits"
        if pipe_conduits_path in hdf_file:
            self._hdf_cache['has_pipe_conduits'] = True
            if 'Attributes' in hdf_file[pipe_conduits_path]:
                self._hdf_cache['pipe_conduits_count'] = len(
                    hdf_file[f"{pipe_conduits_path}/Attributes"][()]
                )
        else:
            self._hdf_cache['has_pipe_conduits'] = False
            
        pipe_nodes_path = "Geometry/Pipe Nodes"
        if pipe_nodes_path in hdf_file:
            self._hdf_cache['has_pipe_nodes'] = True
            if 'Attributes' in hdf_file[pipe_nodes_path]:
                self._hdf_cache['pipe_nodes_count'] = len(
                    hdf_file[f"{pipe_nodes_path}/Attributes"][()]
                )
        else:
            self._hdf_cache['has_pipe_nodes'] = False
        
        # Check for results and get simulation time
        plan_info = hdf_file.get("Plan Data/Plan Information")
        if plan_info and 'Simulation Start Time' in plan_info.attrs:
            time_str = plan_info.attrs['Simulation Start Time']
            self._hdf_cache['simulation_start_time'] = datetime.strptime(
                time_str.decode('utf-8'), "%d%b%Y %H:%M:%S"
            )
            self._hdf_cache['has_results'] = True

    # --- HDF Data Extraction Methods ---

    def _get_breaklines_direct(self, hdf_file, sr):
        """Extracts 2D breaklines from HDF file with optimized numpy operations."""
        try:
            breaklines_path = "Geometry/2D Flow Area Break Lines"
            if breaklines_path not in hdf_file:
                return [], []
            
            bl_line_data = hdf_file[breaklines_path]
            attributes = bl_line_data["Attributes"][()]
            polyline_info = bl_line_data["Polyline Info"][()]
            polyline_points = bl_line_data["Polyline Points"][()]
            
            # Vectorized filtering of valid breaklines
            valid_mask = polyline_info[:, 1] >= 2  # pnt_cnt >= 2
            valid_indices = np.where(valid_mask)[0]
            
            valid_data, geometries = [], []
            
            for idx in valid_indices:
                pnt_start, pnt_cnt, part_start, part_cnt = polyline_info[idx]
                attr_row = attributes[idx]
                
                name = attr_row["Name"]
                name = name.decode('utf-8', 'ignore').strip() if isinstance(name, bytes) else str(name)
                
                try:
                    # Extract points efficiently
                    points = polyline_points[pnt_start:pnt_start + pnt_cnt]
                    
                    if part_cnt == 1:
                        # Single part - direct creation
                        arcpy_array = arcpy.Array([arcpy.Point(p[0], p[1]) for p in points])
                        geom = arcpy.Polyline(arcpy_array, sr)
                    else:
                        # Multi-part polyline
                        parts = bl_line_data["Polyline Parts"][()][part_start:part_start + part_cnt]
                        all_parts_array = arcpy.Array()
                        
                        for part_pnt_start, part_pnt_cnt in parts:
                            if part_pnt_cnt > 1:
                                part_points = points[part_pnt_start:part_pnt_start + part_pnt_cnt]
                                part_array = arcpy.Array([arcpy.Point(p[0], p[1]) for p in part_points])
                                all_parts_array.add(part_array)
                        
                        if all_parts_array.count == 0:
                            continue
                        geom = arcpy.Polyline(all_parts_array, sr)
                    
                    valid_data.append({
                        'bl_id': int(idx),
                        'Name': name,
                        'CellSpaceNear': float(attr_row["Cell Spacing Near"]),
                        'CellSpaceFar': float(attr_row["Cell Spacing Far"]),
                        'NearRepeats': int(attr_row["Near Repeats"]),
                        'ProtectRadius': int(attr_row["Protection Radius"])
                    })
                    geometries.append(geom)
                    
                except Exception as e:
                    arcpy.AddWarning(f"Error processing breakline {idx}: {str(e)}")
                    continue
            
            return valid_data, geometries
            
        except Exception as e:
            arcpy.AddError(f"HDF Read Error (Breaklines): {e}")
            raise arcpy.ExecuteError("Failed to read breaklines from HDF file")

    def _get_bc_lines_direct(self, hdf_file, sr):
        """Extracts 2D boundary condition lines from HDF file."""
        try:
            bc_lines_path = "Geometry/Boundary Condition Lines"
            if bc_lines_path not in hdf_file:
                return [], []
            
            # Get boundary condition line data
            bc_attrs = hdf_file[f"{bc_lines_path}/Attributes"][()]
            polyline_info = hdf_file[f"{bc_lines_path}/Polyline Info"][()]
            polyline_points = hdf_file[f"{bc_lines_path}/Polyline Points"][()]
            
            # Check if multi-part data exists
            has_parts = f"{bc_lines_path}/Polyline Parts" in hdf_file
            if has_parts:
                polyline_parts = hdf_file[f"{bc_lines_path}/Polyline Parts"][()]
            
            # Vectorized filtering of valid boundary condition lines
            valid_mask = polyline_info[:, 1] >= 2  # pnt_cnt >= 2
            valid_indices = np.where(valid_mask)[0]
            
            valid_data, geometries = [], []
            
            for idx in valid_indices:
                pnt_start, pnt_cnt, part_start, part_cnt = polyline_info[idx]
                attr_row = bc_attrs[idx]
                
                # Extract attributes
                name = attr_row["Name"]
                name = name.decode('utf-8', 'ignore').strip() if isinstance(name, bytes) else str(name)
                
                bc_type = attr_row["Type"]
                bc_type = bc_type.decode('utf-8', 'ignore').strip() if isinstance(bc_type, bytes) else str(bc_type)
                
                try:
                    # Extract points efficiently
                    points = polyline_points[pnt_start:pnt_start + pnt_cnt]
                    
                    if part_cnt == 1 or not has_parts:
                        # Single part - direct creation
                        arcpy_array = arcpy.Array([arcpy.Point(p[0], p[1]) for p in points])
                        geom = arcpy.Polyline(arcpy_array, sr)
                    else:
                        # Multi-part polyline
                        parts = polyline_parts[part_start:part_start + part_cnt]
                        all_parts_array = arcpy.Array()
                        
                        for part_pnt_start, part_pnt_cnt in parts:
                            if part_pnt_cnt > 1:
                                part_points = points[part_pnt_start:part_pnt_start + part_pnt_cnt]
                                part_array = arcpy.Array([arcpy.Point(p[0], p[1]) for p in part_points])
                                all_parts_array.add(part_array)
                        
                        if all_parts_array.count == 0:
                            continue
                        geom = arcpy.Polyline(all_parts_array, sr)
                    
                    valid_data.append({
                        'bc_id': int(idx),
                        'Name': name,
                        'Type': bc_type
                    })
                    geometries.append(geom)
                    
                except Exception as e:
                    arcpy.AddWarning(f"Error processing boundary condition line {idx}: {str(e)}")
                    continue
            
            return valid_data, geometries
            
        except Exception as e:
            arcpy.AddError(f"HDF Read Error (Boundary Condition Lines): {e}")
            raise arcpy.ExecuteError("Failed to read boundary condition lines from HDF file")

    def _get_pipe_conduits_direct(self, hdf_file, sr):
        """Extracts pipe conduits from HDF file with dynamic attributes."""
        try:
            conduits_path = "Geometry/Pipe Conduits"
            if conduits_path not in hdf_file:
                return [], []
            
            conduits_group = hdf_file[conduits_path]
            
            # Get attributes
            if 'Attributes' not in conduits_group:
                return [], []
            
            attributes = conduits_group['Attributes'][()]
            
            # Get polyline geometry data
            if 'Polyline Info' not in conduits_group or 'Polyline Points' not in conduits_group:
                return [], []
            
            polyline_info = conduits_group['Polyline Info'][()]
            polyline_points = conduits_group['Polyline Points'][()]
            
            valid_data, geometries = [], []
            
            # Debug: Show original field names
            if len(attributes) > 0:
                arcpy.AddMessage(f"DEBUG: Original HDF field names: {list(attributes.dtype.names)}")
            
            # Process each conduit
            for idx, (info, attr_row) in enumerate(zip(polyline_info, attributes)):
                point_start_idx, point_count = info[0], info[1]
                
                if point_count < 2:
                    continue
                
                try:
                    # Extract points for this conduit
                    coords = polyline_points[point_start_idx:point_start_idx + point_count]
                    
                    # Create polyline geometry
                    arcpy_array = arcpy.Array([arcpy.Point(p[0], p[1]) for p in coords])
                    geom = arcpy.Polyline(arcpy_array, sr)
                    
                    # Build attribute dictionary dynamically
                    attr_dict = {'conduit_id': int(idx)}
                    
                    # Process all attribute fields
                    attr_names = attributes.dtype.names
                    for field_name in attr_names:
                        value = attr_row[field_name]
                        
                        # Decode bytes to string if necessary
                        if isinstance(value, (bytes, np.bytes_)):
                            value = value.decode('utf-8', 'ignore').strip()
                        elif isinstance(value, np.ndarray) and value.dtype.kind == 'S':
                            value = value.tobytes().decode('utf-8', 'ignore').strip()
                        
                        # Clean field name for ArcGIS compatibility
                        clean_name = field_name.replace(' ', '_').replace(':', '_').replace(';', '_').replace(',', '_').replace('(', '_').replace(')', '_').replace("'", '')
                        
                        # Fix known typos in HDF field names
                        if 'Condtui_Connections' in clean_name:
                            clean_name = clean_name.replace('Condtui_Connections', 'Conduit_Connections')
                            if idx == 0:  # Only log for first record
                                arcpy.AddMessage(f"DEBUG: Fixed typo in field name 'Condtui_Connections' to 'Conduit_Connections'")
                        
                        # Special handling for exact "Shape" field (case insensitive)
                        if field_name.upper() == 'SHAPE':
                            clean_name = 'Shape_Type'
                            if idx == 0:
                                arcpy.AddMessage(f"DEBUG: Renamed 'Shape' field to 'Shape_Type' to avoid system field conflict")
                        # Rename other fields that conflict with system fields
                        elif clean_name.upper() in ['OBJECTID', 'SHAPE', 'SHAPE_LENGTH', 'SHAPE_AREA', 'SHAPE_LENG']:
                            original_clean = clean_name
                            clean_name = f"{clean_name}_USER"
                            if idx == 0:  # Only log for first record to avoid spam
                                arcpy.AddMessage(f"DEBUG: Renamed system field '{original_clean}' to '{clean_name}'")
                            
                        attr_dict[clean_name] = value
                    
                    valid_data.append(attr_dict)
                    geometries.append(geom)
                    
                except Exception as e:
                    arcpy.AddWarning(f"Error processing pipe conduit {idx}: {str(e)}")
                    continue
            
            # Debug: Show cleaned field names from first record
            if valid_data:
                arcpy.AddMessage(f"DEBUG: Cleaned field names: {list(valid_data[0].keys())}")
            
            return valid_data, geometries
            
        except Exception as e:
            arcpy.AddError(f"HDF Read Error (Pipe Conduits): {e}")
            raise arcpy.ExecuteError("Failed to read pipe conduits from HDF file")

    def _get_pipe_nodes_direct(self, hdf_file, sr):
        """Extracts pipe nodes from HDF file with dynamic attributes."""
        try:
            nodes_path = "Geometry/Pipe Nodes"
            if nodes_path not in hdf_file:
                return [], []
            
            nodes_group = hdf_file[nodes_path]
            
            # Get attributes
            if 'Attributes' not in nodes_group:
                return [], []
            
            attributes = nodes_group['Attributes'][()]
            
            # Get points data
            if 'Points' not in nodes_group:
                return [], []
            
            points = nodes_group['Points'][()]
            
            valid_data, geometries = [], []
            
            # Debug: Show original field names
            if len(attributes) > 0:
                arcpy.AddMessage(f"DEBUG: Original HDF field names: {list(attributes.dtype.names)}")
            
            # Process each node
            for idx, (xy, attr_row) in enumerate(zip(points, attributes)):
                if len(xy) < 2:
                    continue
                
                try:
                    # Create point geometry
                    geom = arcpy.PointGeometry(arcpy.Point(xy[0], xy[1]), sr)
                    
                    # Build attribute dictionary dynamically
                    attr_dict = {'node_id': int(idx)}
                    
                    # Process all attribute fields
                    attr_names = attributes.dtype.names
                    for field_name in attr_names:
                        value = attr_row[field_name]
                        
                        # Decode bytes to string if necessary
                        if isinstance(value, (bytes, np.bytes_)):
                            value = value.decode('utf-8', 'ignore').strip()
                        elif isinstance(value, np.ndarray) and value.dtype.kind == 'S':
                            value = value.tobytes().decode('utf-8', 'ignore').strip()
                        
                        # Clean field name for ArcGIS compatibility
                        clean_name = field_name.replace(' ', '_').replace(':', '_').replace(';', '_').replace(',', '_').replace('(', '_').replace(')', '_').replace("'", '')
                        
                        # Fix known typos in HDF field names
                        if 'Condtui_Connections' in clean_name:
                            clean_name = clean_name.replace('Condtui_Connections', 'Conduit_Connections')
                            if idx == 0:  # Only log for first record
                                arcpy.AddMessage(f"DEBUG: Fixed typo in field name 'Condtui_Connections' to 'Conduit_Connections'")
                        
                        # Special handling for exact "Shape" field (case insensitive)
                        if field_name.upper() == 'SHAPE':
                            clean_name = 'Shape_Type'
                            if idx == 0:
                                arcpy.AddMessage(f"DEBUG: Renamed 'Shape' field to 'Shape_Type' to avoid system field conflict")
                        # Rename other fields that conflict with system fields
                        elif clean_name.upper() in ['OBJECTID', 'SHAPE', 'SHAPE_LENGTH', 'SHAPE_AREA', 'SHAPE_LENG']:
                            original_clean = clean_name
                            clean_name = f"{clean_name}_USER"
                            if idx == 0:  # Only log for first record to avoid spam
                                arcpy.AddMessage(f"DEBUG: Renamed system field '{original_clean}' to '{clean_name}'")
                            
                        attr_dict[clean_name] = value
                    
                    valid_data.append(attr_dict)
                    geometries.append(geom)
                    
                except Exception as e:
                    arcpy.AddWarning(f"Error processing pipe node {idx}: {str(e)}")
                    continue
            
            # Debug: Show cleaned field names from first record
            if valid_data:
                arcpy.AddMessage(f"DEBUG: Cleaned field names: {list(valid_data[0].keys())}")
            
            return valid_data, geometries
            
        except Exception as e:
            arcpy.AddError(f"HDF Read Error (Pipe Nodes): {e}")
            raise arcpy.ExecuteError("Failed to read pipe nodes from HDF file")

    def _get_dynamic_fields_from_data(self, data_list):
        """Determines field definitions from dynamic attribute data."""
        if not data_list:
            return []
        
        # Get all unique field names and their types
        field_info = {}
        field_lengths = {}  # Track max length for text fields
        
        for record in data_list:
            for field_name, value in record.items():
                # Skip None or NaN values for type detection
                if value is None or (isinstance(value, (float, np.floating)) and np.isnan(value)):
                    continue
                    
                if field_name not in field_info:
                    # Determine field type based on value
                    if isinstance(value, (bool, np.bool_)):
                        field_info[field_name] = "SHORT"  # Use SHORT for boolean
                    elif isinstance(value, (int, np.integer)):
                        field_info[field_name] = "LONG"
                    elif isinstance(value, (float, np.floating)):
                        field_info[field_name] = "DOUBLE"
                    else:
                        field_info[field_name] = "TEXT"
                        field_lengths[field_name] = 0
                
                # Track max length for text fields
                if field_info.get(field_name) == "TEXT" and value is not None:
                    str_value = str(value)
                    field_lengths[field_name] = max(field_lengths.get(field_name, 0), len(str_value))
        
        # Check for fields that might have been skipped due to all NaN values
        all_field_names = set()
        for record in data_list:
            all_field_names.update(record.keys())
        
        # Add any missing fields with default type
        for field_name in all_field_names:
            if field_name not in field_info:
                # Default to DOUBLE for numeric-sounding fields, TEXT otherwise
                if any(keyword in field_name.lower() for keyword in ['elevation', 'area', 'length', 'coefficient', 'offset']):
                    field_info[field_name] = "DOUBLE"
                else:
                    field_info[field_name] = "TEXT"
                    field_lengths[field_name] = 50
        
        # Convert to list of tuples for field creation
        fields = []
        for name, ftype in field_info.items():
            if ftype == "TEXT":
                # Ensure minimum field length of 50, max 255
                length = min(max(field_lengths.get(name, 50), 50), 255)
                fields.append((name, ftype, length))
            else:
                fields.append((name, ftype))
        
        return fields

    def _get_mesh_areas_direct(self, hdf_file, sr):
        """Extracts mesh area perimeters from HDF file."""
        try:
            if not self._hdf_cache['mesh_names']:
                return [], []
            
            raw_data = [{'mesh_name': name} for name in self._hdf_cache['mesh_names']]
            geometries = []
            
            flow_areas_path = "Geometry/2D Flow Areas"
            for mesh_name in self._hdf_cache['mesh_names']:
                perimeter_path = f"{flow_areas_path}/{mesh_name}/Perimeter"
                if perimeter_path in hdf_file:
                    coords = hdf_file[perimeter_path][()]
                    # Create polygon directly from numpy array
                    arcpy_array = arcpy.Array([arcpy.Point(p[0], p[1]) for p in coords])
                    geometries.append(arcpy.Polygon(arcpy_array, sr))
                else:
                    geometries.append(None)
            
            return raw_data, geometries
            
        except Exception as e:
            arcpy.AddError(f"HDF Read Error (Perimeters): {e}")
            raise arcpy.ExecuteError()

    def _get_mesh_cell_points_direct(self, hdf_file, sr):
        """Extracts mesh cell centers using vectorized operations."""
        try:
            if not self._hdf_cache['mesh_names']:
                return [], []
            
            raw_data, geometries = [], []
            
            for mesh_name in self._hdf_cache['mesh_names']:
                cell_centers_path = f"Geometry/2D Flow Areas/{mesh_name}/Cells Center Coordinate"
                if cell_centers_path not in hdf_file:
                    arcpy.AddWarning(f"No cell center data found for mesh '{mesh_name}'")
                    continue
                
                # Read all cell centers at once
                cell_centers = hdf_file[cell_centers_path][()]
                num_cells = len(cell_centers)
                
                # Vectorized data creation
                mesh_data = [{'mesh_name': mesh_name, 'cell_id': i} for i in range(num_cells)]
                raw_data.extend(mesh_data)
                
                # Batch create point geometries
                mesh_geometries = [arcpy.PointGeometry(arcpy.Point(coords[0], coords[1]), sr) 
                                 for coords in cell_centers]
                geometries.extend(mesh_geometries)
            
            return raw_data, geometries
            
        except Exception as e:
            arcpy.AddError(f"HDF Read Error (Cell Points): {e}")
            raise arcpy.ExecuteError()

    def _get_mesh_cell_faces_direct(self, hdf_file, sr):
        """Extracts mesh cell faces with optimized coordinate handling."""
        try:
            if not self._hdf_cache['mesh_names']:
                return [], []
            
            raw_data, geometries = [], []
            
            for mesh_name in self._hdf_cache['mesh_names']:
                try:
                    base = f"Geometry/2D Flow Areas/{mesh_name}"
                    
                    # Load all data at once
                    facepoints_index = hdf_file[f"{base}/Faces FacePoint Indexes"][()]
                    facepoints_coords = hdf_file[f"{base}/FacePoints Coordinate"][()]
                    faces_perim_info = hdf_file[f"{base}/Faces Perimeter Info"][()]
                    faces_perim_values = hdf_file[f"{base}/Faces Perimeter Values"][()]
                    
                    # Process faces in batches
                    for face_id, ((p_a, p_b), (s_row, count)) in enumerate(
                        zip(facepoints_index, faces_perim_info)):
                        
                        # Build coordinate array efficiently
                        if count > 0:
                            coords = np.vstack([
                                facepoints_coords[p_a:p_a+1],
                                faces_perim_values[s_row:s_row + count],
                                facepoints_coords[p_b:p_b+1]
                            ])
                        else:
                            coords = np.vstack([
                                facepoints_coords[p_a:p_a+1],
                                facepoints_coords[p_b:p_b+1]
                            ])
                        
                        # Create polyline
                        arcpy_array = arcpy.Array([arcpy.Point(p[0], p[1]) for p in coords])
                        geometries.append(arcpy.Polyline(arcpy_array, sr))
                        raw_data.append({'mesh_name': mesh_name, 'face_id': face_id})
                    
                except KeyError:
                    arcpy.AddWarning(f"No face data for mesh '{mesh_name}'.")
            
            return raw_data, geometries
            
        except Exception as e:
            arcpy.AddError(f"HDF Read Error (Cell Faces): {e}")
            raise arcpy.ExecuteError()

    def _get_mesh_cells_direct(self, hdf_file, sr, precomputed_faces, messages):
        """
        Optimized mesh cell extraction using numpy arrays and pre-computed lookups.
        """
        try:
            messages.addMessage("Starting optimized cell polygon creation...")
            
            if not self._hdf_cache['mesh_names']:
                return [], []
            
            # Build optimized face lookup with numpy arrays
            face_lookup = {}
            face_arrays = {}  # Store face coordinates as numpy arrays
            
            for i, (face_attr, face_geom) in enumerate(zip(precomputed_faces[0], precomputed_faces[1])):
                mesh_name = face_attr['mesh_name']
                face_id = face_attr['face_id']
                
                if mesh_name not in face_lookup:
                    face_lookup[mesh_name] = {}
                    face_arrays[mesh_name] = {}
                
                face_lookup[mesh_name][face_id] = face_geom
                
                # Store coordinates as numpy array for faster access
                if face_geom and hasattr(face_geom, 'getPart'):
                    part = face_geom.getPart(0)
                    coords = np.array([[part.getObject(i).X, part.getObject(i).Y]
                                     for i in range(part.count) if part.getObject(i)])
                    if len(coords) > 0:
                        face_arrays[mesh_name][face_id] = coords
            
            raw_data, geometries = [], []
            total_cells_processed = 0
            
            for mesh_name in self._hdf_cache['mesh_names']:
                messages.addMessage(f"\nProcessing mesh '{mesh_name}'...")
                
                try:
                    base = f"Geometry/2D Flow Areas/{mesh_name}"
                    
                    # Load cell-face relationships
                    cell_face_info = hdf_file[f"{base}/Cells Face and Orientation Info"][()]
                    cell_face_values = hdf_file[f"{base}/Cells Face and Orientation Values"][()]
                    
                    # Extract as numpy arrays for efficiency
                    face_indices = cell_face_values[:, 0].astype(np.int32)
                    orientations = cell_face_values[:, 1].astype(np.int32)
                    
                    mesh_faces = face_lookup.get(mesh_name, {})
                    mesh_face_arrays = face_arrays.get(mesh_name, {})
                    
                    num_cells = len(cell_face_info)
                    cells_created = 0
                    
                    # Process in batches for progress reporting
                    batch_size = 5000
                    
                    for batch_start in range(0, num_cells, batch_size):
                        batch_end = min(batch_start + batch_size, num_cells)
                        batch_created = 0
                        
                        for cell_id in range(batch_start, batch_end):
                            start, length = cell_face_info[cell_id]
                            
                            if length < 3:  # Need at least 3 faces
                                continue
                            
                            # Get face indices for this cell
                            cell_face_ids = face_indices[start:start + length]
                            cell_orientations = orientations[start:start + length]
                            
                            # Collect face geometries
                            face_geoms = []
                            
                            for j, (face_id, orientation) in enumerate(zip(cell_face_ids, cell_orientations)):
                                if face_id in mesh_faces:
                                    face_geom = mesh_faces[face_id]
                                    
                                    # Handle orientation
                                    if orientation < 0 and face_id in mesh_face_arrays:
                                        # Create reversed geometry using numpy array
                                        coords = mesh_face_arrays[face_id][::-1]
                                        arcpy_array = arcpy.Array([arcpy.Point(x, y) for x, y in coords])
                                        face_geom = arcpy.Polyline(arcpy_array, sr)
                                    
                                    if face_geom:
                                        face_geoms.append(face_geom)
                            
                            if len(face_geoms) >= 3:
                                # Use optimized polygon construction
                                polygon = polygonize_arcpy_optimized(face_geoms, sr)
                                
                                if polygon and polygon.area > 0:
                                    raw_data.append({'mesh_name': mesh_name, 'cell_id': cell_id})
                                    geometries.append(polygon)
                                    cells_created += 1
                                    batch_created += 1
                                    total_cells_processed += 1
                        
                        # Progress update
                        if batch_end % 10000 == 0 or batch_end == num_cells:
                            messages.addMessage(
                                f"  Processed {batch_end}/{num_cells} cells in mesh '{mesh_name}' "
                                f"({cells_created} valid polygons)"
                            )
                    
                    messages.addMessage(
                        f"  Completed mesh '{mesh_name}': {cells_created} cells created"
                    )
                    
                except Exception as e:
                    messages.addErrorMessage(f"Error processing mesh '{mesh_name}': {str(e)}")
                    continue
            
            messages.addMessage(f"\nTotal cells processed: {total_cells_processed}")
            return raw_data, geometries
            
        except Exception as e:
            messages.addErrorMessage(f"Fatal error in cell creation: {str(e)}")
            raise arcpy.ExecuteError()

    def _get_max_wse_points_direct(self, hdf_file, sr):
        """Extracts maximum water surface elevation points with vectorized operations."""
        try:
            if not self._hdf_cache['has_results']:
                arcpy.AddError("No results data found in HDF file.")
                return [], []
            
            raw_data, geometries = [], []
            start_time = self._hdf_cache['simulation_start_time']
            
            arcpy.AddMessage(f'Simulation start time: {start_time}')
            
            for mesh_name in self._hdf_cache['mesh_names']:
                arcpy.AddMessage(f'Processing max WSE for mesh: {mesh_name}')
                
                # Get cell centers
                centers_path = f"Geometry/2D Flow Areas/{mesh_name}/Cells Center Coordinate"
                if centers_path not in hdf_file:
                    arcpy.AddWarning(f"No cell centers found for mesh '{mesh_name}'. Skipping.")
                    continue
                
                cell_centers = hdf_file[centers_path][()]
                
                # Get maximum water surface data
                summary_path = f"Results/Unsteady/Output/Output Blocks/Base Output/Summary Output/2D Flow Areas/{mesh_name}/Maximum Water Surface"
                if summary_path not in hdf_file:
                    arcpy.AddWarning(f"No 'Maximum Water Surface' data for mesh '{mesh_name}'. Skipping.")
                    continue
                
                max_wse_data = hdf_file[summary_path][:]
                
                # Data is 2D array: row 0 = values, row 1 = times (in days)
                if max_wse_data.ndim == 2 and max_wse_data.shape[0] == 2:
                    wse_values = max_wse_data[0, :]
                    time_in_days = max_wse_data[1, :]
                else:
                    arcpy.AddWarning(f"Unexpected data format for mesh '{mesh_name}'. Skipping.")
                    continue
                
                num_cells = min(len(cell_centers), len(wse_values))
                
                # Vectorized time conversion
                time_deltas = np.array([timedelta(days=float(t)) for t in time_in_days[:num_cells]])
                times_of_max = [start_time + td for td in time_deltas]
                
                # Create data and geometries in batches
                mesh_data = [{
                    'mesh_name': mesh_name,
                    'cell_id': i,
                    'max_wse': float(wse_values[i]),
                    'max_wse_time': times_of_max[i]
                } for i in range(num_cells)]
                
                mesh_geometries = [
                    arcpy.PointGeometry(arcpy.Point(cell_centers[i][0], cell_centers[i][1]), sr)
                    for i in range(num_cells)
                ]
                
                raw_data.extend(mesh_data)
                geometries.extend(mesh_geometries)
            
            return raw_data, geometries
            
        except Exception as e:
            arcpy.AddError(f"HDF Read Error (Max WSE Points): {e}")
            raise arcpy.ExecuteError()

    def _get_max_face_velocity_points_direct(self, hdf_file, sr):
        """Extracts maximum face velocity points with optimized centroid calculation."""
        try:
            if not self._hdf_cache['has_results']:
                arcpy.AddError("No results data found in HDF file.")
                return [], []
            
            raw_data, geometries = [], []
            start_time = self._hdf_cache['simulation_start_time']
            
            arcpy.AddMessage(f'Simulation start time: {start_time}')
            
            # Get face geometries
            face_data, face_geoms = self._get_mesh_cell_faces_direct(hdf_file, sr)
            
            # Build lookup
            face_lookup = {}
            for i, (face_attr, face_geom) in enumerate(zip(face_data, face_geoms)):
                mesh_name = face_attr['mesh_name']
                face_id = face_attr['face_id']
                
                if mesh_name not in face_lookup:
                    face_lookup[mesh_name] = {}
                
                face_lookup[mesh_name][face_id] = face_geom
            
            for mesh_name in self._hdf_cache['mesh_names']:
                arcpy.AddMessage(f'Processing max face velocity for mesh: {mesh_name}')
                
                # Get maximum face velocity data
                summary_path = f"Results/Unsteady/Output/Output Blocks/Base Output/Summary Output/2D Flow Areas/{mesh_name}/Maximum Face Velocity"
                if summary_path not in hdf_file:
                    arcpy.AddWarning(f"No 'Maximum Face Velocity' data for mesh '{mesh_name}'. Skipping.")
                    continue
                
                max_vel_data = hdf_file[summary_path][:]
                
                # Data is 2D array: row 0 = values, row 1 = times (in days)
                if max_vel_data.ndim == 2 and max_vel_data.shape[0] == 2:
                    vel_values = max_vel_data[0, :]
                    time_in_days = max_vel_data[1, :]
                else:
                    arcpy.AddWarning(f"Unexpected data format for mesh '{mesh_name}'. Skipping.")
                    continue
                
                # Process faces
                mesh_faces = face_lookup.get(mesh_name, {})
                
                for face_id in range(len(vel_values)):
                    if face_id not in mesh_faces:
                        continue
                    
                    face_geom = mesh_faces[face_id]
                    if not face_geom:
                        continue
                    
                    # Calculate centroid using optimized function
                    centroid_pt = get_polyline_centroid_vectorized(face_geom)
                    if not centroid_pt:
                        continue
                    
                    # Convert time
                    time_of_max = start_time + timedelta(days=float(time_in_days[face_id]))
                    
                    raw_data.append({
                        'mesh_name': mesh_name,
                        'face_id': face_id,
                        'max_vel': float(vel_values[face_id]),
                        'time_of_max': time_of_max
                    })
                    
                    geometries.append(arcpy.PointGeometry(centroid_pt, sr))
            
            return raw_data, geometries
            
        except Exception as e:
            arcpy.AddError(f"HDF Read Error (Max Face Velocity Points): {e}")
            raise arcpy.ExecuteError()

    # --- Optimized Feature Writing Method ---

    def _write_features_to_fc(self, output_fc, sr, geom_type, fields, data, geometries, messages):
        """Optimized feature writing with batch operations."""
        total_features = len(geometries)
        if total_features == 0:
            messages.addWarningMessage(f"No features found for {os.path.basename(output_fc)}. Layer will not be created.")
            return

        messages.addMessage(f"Creating feature class: {os.path.basename(output_fc)} ({total_features} features)")
        
        try:
            output_path, output_name = os.path.split(output_fc)
            if arcpy.Exists(output_fc):
                arcpy.management.Delete(output_fc)
            arcpy.management.CreateFeatureclass(output_path, output_name, geom_type, spatial_reference=sr)
            
            field_names = []
            for field_def in fields:
                field_name = field_def[0]
                field_type = field_def[1]
                
                if field_type == "DATE":
                    arcpy.management.AddField(output_fc, field_name, "DATE")
                elif field_type == "TEXT" and len(field_def) > 2:
                    # Use the calculated field length for text fields
                    field_length = field_def[2]
                    arcpy.management.AddField(output_fc, field_name, field_type, field_length=field_length)
                else:
                    arcpy.management.AddField(output_fc, field_name, field_type)
                field_names.append(field_name)
                
        except Exception as e:
            messages.addErrorMessage(f"Failed to create output feature class {output_fc}: {e}")
            return
        
        try:
            field_names_with_shape = ["SHAPE@"] + field_names
            features_inserted = 0
            
            # Prepare all rows at once
            all_rows = []
            missing_fields_warned = set()  # Track which fields we've already warned about
            
            for i, geom in enumerate(geometries):
                if geom is None or geom.pointCount == 0:
                    continue
                
                row_data = data[i]
                row_values = []
                for name in field_names:
                    # Use the exact field name from the data, accounting for any cleanup
                    value = row_data.get(name)
                    if value is None and name not in row_data:
                        # Only warn once per field across all rows
                        if name not in missing_fields_warned:
                            messages.addWarning(f"Field '{name}' not found in data. Available fields: {list(row_data.keys())}")
                            missing_fields_warned.add(name)
                        value = None  # Use None for missing fields
                    
                    # Handle special numpy types and NaN values
                    if value is not None:
                        # Convert numpy NaN to None for proper NULL handling
                        if isinstance(value, (float, np.floating)) and np.isnan(value):
                            value = None
                        # Convert numpy types to Python types
                        elif isinstance(value, np.integer):
                            value = int(value)
                        elif isinstance(value, np.floating):
                            value = float(value)
                        elif isinstance(value, np.bool_):
                            value = int(value)  # Convert bool to 0/1 for SHORT field
                    
                    row_values.append(value)
                all_rows.append([geom] + row_values)
            
            # Batch insert with larger chunks
            batch_size = 10000
            with arcpy.da.InsertCursor(output_fc, field_names_with_shape) as cursor:
                for i in range(0, len(all_rows), batch_size):
                    batch = all_rows[i:i + batch_size]
                    for row in batch:
                        cursor.insertRow(row)
                    features_inserted += len(batch)
                    
                    if features_inserted % 50000 == 0:
                        messages.addMessage(f"  Inserted {features_inserted}/{len(all_rows)} features...")
            
            messages.addMessage(f"Successfully created {os.path.basename(output_fc)} with {features_inserted} features.")
            
        except Exception as e:
            # If the error is about a missing field, show more diagnostic info
            if "Cannot find field" in str(e) or "not found in data" in str(e):
                messages.addErrorMessage(f"Field mapping error for {output_fc}:")
                messages.addErrorMessage(f"  Expected fields: {field_names}")
                if all_rows and len(all_rows) > 0:
                    sample_data = data[0] if data else {}
                    messages.addErrorMessage(f"  Available fields in data: {list(sample_data.keys())}")
                messages.addErrorMessage(f"  Original error: {e}")
            else:
                messages.addErrorMessage(f"Failed during feature creation for {output_fc}: {e}")

    # --- Main Execution Logic ---
    def execute(self, parameters, messages):
        hdf_path = parameters[0].valueAsText
        
        # Get selected elements
        geometry_elements = parameters[2].values if parameters[2].values else []
        results_elements = parameters[3].values if parameters[3].values else []
        elements_to_load = (geometry_elements or []) + (results_elements or [])
        
        if not elements_to_load:
            messages.addErrorMessage("No elements selected for loading. Please select at least one element.")
            raise arcpy.ExecuteError
        
        # Get projection
        proj_wkt = get_ras_projection_wkt(hdf_path)
        sr = None
        if proj_wkt:
            sr = arcpy.SpatialReference()
            sr.loadFromString(proj_wkt)
            messages.addMessage(f"CRS '{sr.name}' found in HEC-RAS project files.")
        elif parameters[1].value:
            sr = parameters[1].value
            messages.addMessage(f"Using user-defined override CRS: {sr.name}")
        else:
            messages.addErrorMessage("CRS could not be determined. Please use the Override CRS parameter.")
            raise arcpy.ExecuteError
        
        # Open HDF file once and cache metadata
        with h5py.File(hdf_path, 'r') as hdf_file:
            messages.addMessage("Caching HDF metadata...")
            self._cache_hdf_metadata(hdf_file)
            
            # Pre-compute faces if needed
            precomputed_faces = None
            if (self.CELL_FACES in elements_to_load or 
                self.CELL_POLYS in elements_to_load or 
                self.MAX_FACE_VEL_POINTS in elements_to_load):
                messages.addMessage("Pre-loading Mesh Cell Faces...")
                precomputed_faces = self._get_mesh_cell_faces_direct(hdf_file, sr)
            
            # Process geometry elements
            if self.BREAKLINES in elements_to_load and parameters[4].valueAsText:
                output_fc = parameters[4].valueAsText
                data, geoms = self._get_breaklines_direct(hdf_file, sr)
                fields = [("bl_id", "LONG"), ("Name", "TEXT"), ("CellSpaceNear", "FLOAT"), 
                         ("CellSpaceFar", "FLOAT"), ("NearRepeats", "LONG"), ("ProtectRadius", "LONG")]
                self._write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
            
            if self.BC_LINES in elements_to_load and parameters[5].valueAsText:
                output_fc = parameters[5].valueAsText
                data, geoms = self._get_bc_lines_direct(hdf_file, sr)
                fields = [("bc_id", "LONG"), ("Name", "TEXT"), ("Type", "TEXT")]
                self._write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
            
            if self.PERIMETERS in elements_to_load and parameters[6].valueAsText:
                output_fc = parameters[6].valueAsText
                data, geoms = self._get_mesh_areas_direct(hdf_file, sr)
                fields = [("mesh_name", "TEXT")]
                self._write_features_to_fc(output_fc, sr, "POLYGON", fields, data, geoms, messages)
            
            if self.CELL_POINTS in elements_to_load and parameters[7].valueAsText:
                output_fc = parameters[7].valueAsText
                data, geoms = self._get_mesh_cell_points_direct(hdf_file, sr)
                fields = [("mesh_name", "TEXT"), ("cell_id", "LONG")]
                self._write_features_to_fc(output_fc, sr, "POINT", fields, data, geoms, messages)
            
            if self.CELL_FACES in elements_to_load and parameters[8].valueAsText:
                output_fc = parameters[8].valueAsText
                if precomputed_faces:
                    data, geoms = precomputed_faces
                    fields = [("mesh_name", "TEXT"), ("face_id", "LONG")]
                    self._write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
            
            if self.CELL_POLYS in elements_to_load and parameters[9].valueAsText:
                output_fc = parameters[9].valueAsText
                if precomputed_faces:
                    messages.addMessage("Constructing cell polygons from faces...")
                    data, geoms = self._get_mesh_cells_direct(hdf_file, sr, precomputed_faces, messages)
                    fields = [("mesh_name", "TEXT"), ("cell_id", "LONG")]
                    self._write_features_to_fc(output_fc, sr, "POLYGON", fields, data, geoms, messages)
            
            # Process results elements
            if self.MAX_WSE_POINTS in elements_to_load and parameters[10].valueAsText:
                output_fc = parameters[10].valueAsText
                messages.addMessage("Extracting Maximum Water Surface Elevation points...")
                data, geoms = self._get_max_wse_points_direct(hdf_file, sr)
                fields = [("mesh_name", "TEXT"), ("cell_id", "LONG"), ("max_wse", "DOUBLE"), 
                         ("max_wse_time", "DATE")]
                self._write_features_to_fc(output_fc, sr, "POINT", fields, data, geoms, messages)
            
            if self.MAX_FACE_VEL_POINTS in elements_to_load and parameters[11].valueAsText:
                output_fc = parameters[11].valueAsText
                messages.addMessage("Extracting Maximum Face Velocity points...")
                data, geoms = self._get_max_face_velocity_points_direct(hdf_file, sr)
                fields = [("mesh_name", "TEXT"), ("face_id", "LONG"), ("max_vel", "DOUBLE"), 
                         ("time_of_max", "DATE")]
                self._write_features_to_fc(output_fc, sr, "POINT", fields, data, geoms, messages)
            
            # Process pipe network elements
            if self.PIPE_CONDUITS in elements_to_load and parameters[12].valueAsText:
                output_fc = parameters[12].valueAsText
                messages.addMessage("Extracting Pipe Conduits...")
                data, geoms = self._get_pipe_conduits_direct(hdf_file, sr)
                if data:
                    # Get dynamic fields from the data
                    fields = self._get_dynamic_fields_from_data(data)
                    self._write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
                else:
                    messages.addWarning("No pipe conduits found in the HDF file.")
            
            if self.PIPE_NODES in elements_to_load and parameters[13].valueAsText:
                output_fc = parameters[13].valueAsText
                messages.addMessage("Extracting Pipe Nodes...")
                data, geoms = self._get_pipe_nodes_direct(hdf_file, sr)
                if data:
                    # Get dynamic fields from the data
                    fields = self._get_dynamic_fields_from_data(data)
                    self._write_features_to_fc(output_fc, sr, "POINT", fields, data, geoms, messages)
                else:
                    messages.addWarning("No pipe nodes found in the HDF file.")
        
        messages.addMessage("\nProcessing complete.")
        return
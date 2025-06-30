# -*- coding: utf-8 -*-
"""
LoadHECRAS2DGeometry.py

Tool for loading HEC-RAS 2D geometry layers from HDF files including mesh elements,
breaklines, boundary conditions, and pipe networks.
"""

import arcpy
import os
import h5py
import numpy as np
from collections import defaultdict

# Import helper functions from utils
from .utils import (
    get_ras_projection_wkt,
    polygonize_arcpy_optimized,
    get_polyline_centroid_vectorized,
    cache_hdf_metadata,
    write_features_to_fc,
    get_dynamic_fields_from_data
)


class LoadHECRAS2DGeometry(object):
    """
    Loads 2D geometry elements from a HEC-RAS HDF file.
    """
    def __init__(self):
        self.label = "Load HEC-RAS 2D Geometry Layers"
        self.description = """Extracts 2D geometry elements from a HEC-RAS HDF file including mesh elements, breaklines, boundary conditions, and pipe networks.
        
        This tool extracts various 2D geometry elements from HEC-RAS geometry (g*.hdf) or plan (p*.hdf) files.
        
        Available geometry elements include:
        • 2D Breaklines - Mesh refinement lines with cell spacing attributes
        • 2D Boundary Condition Lines - External and internal boundary conditions
        • Mesh Area Perimeters - 2D flow area boundaries
        • Mesh Cell Centers - Point locations at the center of each mesh cell
        • Mesh Cell Faces - Line geometries representing cell edges
        • Mesh Cells (Polygons) - Full polygon representation of mesh cells
        • Pipe Conduits - Storm/sewer pipe networks (if present)
        • Pipe Nodes - Junction points in pipe networks (if present)
        
        Note: Mesh cell polygon creation can be time-consuming for large meshes (>100,000 cells)."""
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
        
        # Cache for HDF metadata
        self._hdf_cache = {}

    def getParameterInfo(self):
        geometry_elements = [self.BREAKLINES, self.BC_LINES, self.PERIMETERS, self.CELL_POINTS, 
                           self.CELL_FACES, self.CELL_POLYS, self.PIPE_CONDUITS, self.PIPE_NODES, self.PIPE_NETWORKS]

        params = [
            arcpy.Parameter(displayName="Geometry or Plan HDF File", name="input_hdf", datatype="DEFile", 
                          parameterType="Required", direction="Input"),
            arcpy.Parameter(displayName="Override CRS (Optional)", name="override_crs", datatype="GPSpatialReference", 
                          parameterType="Optional", direction="Input"),
            
            # Geometry elements to load
            arcpy.Parameter(displayName="Geometry Elements to Load", name="geometry_elements", datatype="GPString", 
                          parameterType="Required", direction="Input", multiValue=True),
            
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
            arcpy.Parameter(displayName="Output Pipe Conduits", name="output_pipe_conduits", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output Pipe Nodes", name="output_pipe_nodes", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output Pipe Networks", name="output_pipe_networks", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs")
        ]
        
        # Configure HDF file filter
        params[0].filter.list = ["hdf", "g*.hdf", "p*.hdf"]
        params[0].description = "Select a HEC-RAS geometry file (g*.hdf) or plan file (p*.hdf) containing 2D geometry data."
        
        params[1].description = """Specify a coordinate reference system if it cannot be determined from the HEC-RAS project files. 
        The tool will first attempt to read the CRS from the HDF file or associated .prj files."""
        
        # Set filters for multi-value parameters
        params[2].filter.type = "ValueList"
        params[2].filter.list = geometry_elements
        params[2].value = [self.PERIMETERS]  # Default selection
        params[2].description = """Select one or more geometry elements to extract from the HDF file. 
        Each selected element will create a separate output feature class."""
        
        # Set default output paths and descriptions
        params[3].value = r"memory\Breaklines"
        params[3].description = "Output feature class for 2D breaklines with cell spacing attributes."
        
        params[4].value = r"memory\BoundaryConditionLines"
        params[4].description = "Output feature class for 2D boundary condition lines."
        
        params[5].value = r"memory\MeshPerimeters"
        params[5].description = "Output feature class for 2D flow area perimeter polygons."
        
        params[6].value = r"memory\MeshCellCenters"
        params[6].description = "Output feature class for mesh cell center points."
        
        params[7].value = r"memory\MeshCellFaces"
        params[7].description = "Output feature class for mesh cell face polylines."
        
        params[8].value = r"memory\MeshCellPolygons"
        params[8].description = "Output feature class for mesh cell polygons. Note: Can be slow for large meshes."
        
        params[9].value = r"memory\PipeConduits"
        params[9].description = "Output feature class for pipe conduits (storm/sewer networks)."
        
        params[10].value = r"memory\PipeNodes"
        params[10].description = "Output feature class for pipe junction nodes."
        
        params[11].value = r"memory\PipeNetworks"
        params[11].description = "Output feature class for pipe network elements."
        
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal validation."""
        # Enable/disable output parameters based on selected elements
        if parameters[2].value:
            selected = parameters[2].valueAsText.split(';') if parameters[2].valueAsText else []
            
            # Enable/disable outputs based on selection
            parameters[3].enabled = self.BREAKLINES in selected
            parameters[4].enabled = self.BC_LINES in selected
            parameters[5].enabled = self.PERIMETERS in selected
            parameters[6].enabled = self.CELL_POINTS in selected
            parameters[7].enabled = self.CELL_FACES in selected
            parameters[8].enabled = self.CELL_POLYS in selected
            parameters[9].enabled = self.PIPE_CONDUITS in selected
            parameters[10].enabled = self.PIPE_NODES in selected
            parameters[11].enabled = self.PIPE_NETWORKS in selected
        return
    
    def updateMessages(self, parameters):
        """Modify messages created by internal validation."""
        # Add warning for large mesh polygon creation
        if parameters[2].value and self.CELL_POLYS in str(parameters[2].value):
            parameters[8].setWarningMessage(
                "Creating cell polygons can be time-consuming for large meshes (>100,000 cells). "
                "Consider using cell centers or faces for visualization instead."
            )
        return

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

    # --- Main Execution Logic ---
    def execute(self, parameters, messages):
        hdf_path = parameters[0].valueAsText
        
        # Get selected elements
        geometry_elements = parameters[2].values if parameters[2].values else []
        
        if not geometry_elements:
            messages.addErrorMessage("No geometry elements selected for loading. Please select at least one element.")
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
            self._hdf_cache = cache_hdf_metadata(hdf_file)
            
            # Pre-compute faces if needed
            precomputed_faces = None
            if (self.CELL_FACES in geometry_elements or 
                self.CELL_POLYS in geometry_elements):
                messages.addMessage("Pre-loading Mesh Cell Faces...")
                precomputed_faces = self._get_mesh_cell_faces_direct(hdf_file, sr)
            
            # Process geometry elements
            if self.BREAKLINES in geometry_elements and parameters[3].valueAsText:
                output_fc = parameters[3].valueAsText
                data, geoms = self._get_breaklines_direct(hdf_file, sr)
                fields = [("bl_id", "LONG"), ("Name", "TEXT"), ("CellSpaceNear", "FLOAT"), 
                         ("CellSpaceFar", "FLOAT"), ("NearRepeats", "LONG"), ("ProtectRadius", "LONG")]
                write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
            
            if self.BC_LINES in geometry_elements and parameters[4].valueAsText:
                output_fc = parameters[4].valueAsText
                data, geoms = self._get_bc_lines_direct(hdf_file, sr)
                fields = [("bc_id", "LONG"), ("Name", "TEXT"), ("Type", "TEXT")]
                write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
            
            if self.PERIMETERS in geometry_elements and parameters[5].valueAsText:
                output_fc = parameters[5].valueAsText
                data, geoms = self._get_mesh_areas_direct(hdf_file, sr)
                fields = [("mesh_name", "TEXT")]
                write_features_to_fc(output_fc, sr, "POLYGON", fields, data, geoms, messages)
            
            if self.CELL_POINTS in geometry_elements and parameters[6].valueAsText:
                output_fc = parameters[6].valueAsText
                data, geoms = self._get_mesh_cell_points_direct(hdf_file, sr)
                fields = [("mesh_name", "TEXT"), ("cell_id", "LONG")]
                write_features_to_fc(output_fc, sr, "POINT", fields, data, geoms, messages)
            
            if self.CELL_FACES in geometry_elements and parameters[7].valueAsText:
                output_fc = parameters[7].valueAsText
                if precomputed_faces:
                    data, geoms = precomputed_faces
                    fields = [("mesh_name", "TEXT"), ("face_id", "LONG")]
                    write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
            
            if self.CELL_POLYS in geometry_elements and parameters[8].valueAsText:
                output_fc = parameters[8].valueAsText
                if precomputed_faces:
                    messages.addMessage("Constructing cell polygons from faces...")
                    data, geoms = self._get_mesh_cells_direct(hdf_file, sr, precomputed_faces, messages)
                    fields = [("mesh_name", "TEXT"), ("cell_id", "LONG")]
                    write_features_to_fc(output_fc, sr, "POLYGON", fields, data, geoms, messages)
            
            # Process pipe network elements
            if self.PIPE_CONDUITS in geometry_elements and parameters[9].valueAsText:
                output_fc = parameters[9].valueAsText
                messages.addMessage("Extracting Pipe Conduits...")
                data, geoms = self._get_pipe_conduits_direct(hdf_file, sr)
                if data:
                    # Get dynamic fields from the data
                    fields = get_dynamic_fields_from_data(data)
                    write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
                else:
                    messages.addWarning("No pipe conduits found in the HDF file.")
            
            if self.PIPE_NODES in geometry_elements and parameters[10].valueAsText:
                output_fc = parameters[10].valueAsText
                messages.addMessage("Extracting Pipe Nodes...")
                data, geoms = self._get_pipe_nodes_direct(hdf_file, sr)
                if data:
                    # Get dynamic fields from the data
                    fields = get_dynamic_fields_from_data(data)
                    write_features_to_fc(output_fc, sr, "POINT", fields, data, geoms, messages)
                else:
                    messages.addWarning("No pipe nodes found in the HDF file.")
        
        messages.addMessage("\nProcessing complete.")
        return
    
    def getHelp(self, tool_name):
        """Return help documentation URL for the tool."""
        help_file = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 
                                "Doc", "RASCommander_Help.html")
        if os.path.exists(help_file):
            return f"file:///{help_file.replace(os.sep, '/')}#load-hec-ras-2d-geometry-layers"
        return None
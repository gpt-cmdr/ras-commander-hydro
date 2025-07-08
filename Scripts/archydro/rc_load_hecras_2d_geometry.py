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
from rc_utils import (
    get_ras_projection_wkt,
    polygonize_arcpy_optimized,
    get_polyline_centroid_vectorized,
    cache_hdf_metadata,
    write_features_to_fc,
    get_dynamic_fields_from_data,
    setup_geodatabase_output,
    get_unique_fc_name,
    add_feature_class_metadata,
    extract_project_and_plan_info,
    create_geodatabase_from_hdf,
    get_feature_dataset_name,
    get_feature_class_name
)


class LoadHECRAS2DGeometry(object):
    """
    Loads 2D geometry elements from a HEC-RAS HDF file.
    """
    def __init__(self):
        # Core properties
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
        
        # Extended metadata properties
        self.summary = "Extract 2D mesh geometry and pipe networks from HEC-RAS HDF files"
        self.usage = """Select a HEC-RAS geometry or plan HDF file and choose which 2D geometry elements to extract.
        
        Steps:
        1. Browse to a HEC-RAS geometry (g*.hdf) or plan (p*.hdf) file
        2. Select which 2D geometry elements to extract
        3. Specify output locations for each selected element
        4. Optionally create an organized geodatabase
        
        Performance considerations:
        • Mesh polygon creation can be slow for meshes > 100,000 cells
        • Consider extracting only cell centers/faces for large meshes
        • Use geodatabase output for better performance with large datasets"""
        
        # Tool behavior
        self.canRunInBackground = False
        self.category = "HEC-RAS Data Import"
        
        # Documentation and credits
        self.tags = ["HEC-RAS", "2D Geometry", "Mesh", "Breaklines", "Boundary Conditions", 
                     "Pipe Networks", "Storm Sewer", "Arc Hydro"]
        self.credits = "CLB Engineering Corporation"
        self.author = "CLB Engineering Corporation"
        self.version = "1.0.0"
        
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
                          parameterType="Optional", direction="Output", category="Outputs"),
            
            # Geodatabase organization parameters
            arcpy.Parameter(displayName="Output Geodatabase (Optional)", name="output_gdb", datatype="DEWorkspace", 
                          parameterType="Optional", direction="Output", category="Output"),
            arcpy.Parameter(displayName="Create New Geodatabase", name="create_gdb", datatype="GPBoolean", 
                          parameterType="Optional", direction="Input", category="Output")
        ]
        
        # Configure HDF file filter
        params[0].filter.list = ["hdf", "g*.hdf", "p*.hdf"]
        params[0].description = """Select a HEC-RAS geometry file (g*.hdf) or plan file (p*.hdf) containing 2D geometry data.
        
        The tool will automatically detect available 2D flow areas and pipe networks in the file.
        
        Supported file types:
        • Geometry files (g01.hdf, g02.hdf, etc.)
        • Plan files with geometry (p01.hdf, p02.hdf, etc.)"""
        params[0].category = "Input Data"
        
        params[1].description = """Specify a coordinate reference system if it cannot be determined from the HEC-RAS project files.
        
        The tool will first attempt to read the CRS from:
        1. The HDF file metadata
        2. Associated .prj files in the project directory
        3. The RAS Mapper projection file
        
        Only provide this parameter if automatic detection fails."""
        params[1].category = "Input Data"
        
        # Set filters for multi-value parameters
        params[2].filter.type = "ValueList"
        params[2].filter.list = geometry_elements
        params[2].value = [self.PERIMETERS]  # Default selection
        params[2].description = """Select one or more geometry elements to extract from the HDF file.
        
        Mesh Elements:
        • 2D Breaklines - Enforce mesh refinement along important features
        • 2D Boundary Condition Lines - Define inflow/outflow boundaries
        • Mesh Area Perimeters - 2D flow area boundaries
        • Mesh Cell Centers - Point at center of each computational cell
        • Mesh Cell Faces - Lines representing cell edges
        • Mesh Cells (Polygons) - Full polygon cells (slow for large meshes)
        
        Pipe Networks:
        • Pipe Conduits - Storm/sewer pipe segments
        • Pipe Nodes - Manholes and junctions
        • Pipe Networks - Complete network elements
        
        Each selected element creates a separate output feature class."""
        params[2].category = "Geometry Selection"
        
        # Set default output paths and descriptions
        params[3].value = r"memory\Breaklines"
        params[3].description = """Output feature class for 2D breaklines.
        
        Attributes include:
        • Name and type
        • Cell spacing along breakline
        • 2D flow area association"""
        
        params[4].value = r"memory\BoundaryConditionLines"
        params[4].description = """Output feature class for 2D boundary condition lines.
        
        Attributes include:
        • BC type (Flow, Stage, Normal Depth, etc.)
        • BC name
        • 2D flow area association"""
        
        params[5].value = r"memory\MeshPerimeters"
        params[5].description = """Output feature class for 2D flow area perimeter polygons.
        
        Attributes include:
        • 2D area name
        • Cell count
        • Minimum cell size
        • Area in acres"""
        
        params[6].value = r"memory\MeshCellCenters"
        params[6].description = """Output feature class for mesh cell center points.
        
        Attributes include:
        • Cell ID
        • 2D flow area name
        • Cell area
        • Elevation (if available)"""
        
        params[7].value = r"memory\MeshCellFaces"
        params[7].description = """Output feature class for mesh cell face polylines.
        
        Represents the edges between computational cells.
        Useful for understanding mesh connectivity."""
        
        params[8].value = r"memory\MeshCellPolygons"
        params[8].description = """Output feature class for mesh cell polygons.
        
        WARNING: Polygon creation can be very slow for large meshes (>100,000 cells).
        Consider using cell centers and faces for large models.
        
        Attributes include:
        • Cell ID
        • Cell area
        • 2D flow area name"""
        
        params[9].value = r"memory\PipeConduits"
        params[9].description = """Output feature class for pipe conduits.
        
        Storm/sewer pipe segments with attributes:
        • Pipe name and material
        • Diameter/dimensions
        • Upstream/downstream nodes
        • Length and slope"""
        
        params[10].value = r"memory\PipeNodes"
        params[10].description = """Output feature class for pipe junction nodes.
        
        Manholes and junctions with attributes:
        • Node name
        • Rim elevation
        • Invert elevation
        • Node type"""
        
        params[11].value = r"memory\PipeNetworks"
        params[11].description = """Output feature class for complete pipe network elements.
        
        Combined pipe network geometry including both conduits and nodes."""
        
        # Geodatabase parameters
        params[12].description = """Specify a geodatabase to organize all output feature classes.
        
        If provided:
        • All outputs will be created in feature datasets within this geodatabase
        • Feature datasets will be organized by geometry type (Mesh, Pipes, etc.)
        • Automatic naming conventions will be applied
        
        Leave empty to use default output locations."""
        params[12].category = "Output Organization"
        
        params[13].value = True  # Default to creating new geodatabase
        params[13].description = """Create a new geodatabase based on the HDF file name.
        
        When enabled:
        • Creates geodatabase named: ProjectName.pXX.gdb
        • Organizes outputs in feature datasets by type
        • Maintains HEC-RAS project structure
        • Optimizes performance for large datasets
        
        Recommended for organizing complete 2D models."""
        params[13].category = "Output Organization"
        
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
            
        # Auto-populate geodatabase path when HDF file is selected
        if parameters[0].value and parameters[0].altered:  # input_hdf
            hdf_path = parameters[0].valueAsText
            
            # If create_gdb is True, auto-populate geodatabase path
            if parameters[13].value:  # create_gdb
                project_name, plan_number, base_name = extract_project_and_plan_info(hdf_path)
                gdb_name = f"{base_name}.gdb"
                gdb_path = os.path.join(os.path.dirname(hdf_path), gdb_name)
                parameters[12].value = gdb_path
        return
    
    def updateMessages(self, parameters):
        """Modify messages created by internal validation."""
        # Add warning for large mesh polygon creation
        if parameters[2].value and self.CELL_POLYS in str(parameters[2].value):
            parameters[8].setWarningMessage(
                "Creating cell polygons can be time-consuming for large meshes (>100,000 cells). "
                "Consider using cell centers or faces for visualization instead."
            )
        
        # Clear geodatabase validation error if create_gdb is True
        if parameters[13].value and parameters[12].hasError():  # create_gdb and output_gdb has error
            parameters[12].clearMessage()
        
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
            arcpy.AddWarning(f"Could not read pipe conduits: {e}")
            return [], []

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
            arcpy.AddWarning(f"Could not read pipe nodes: {e}")
            return [], []

    def _get_pipe_networks_direct(self, hdf_file, sr):
        """Extracts pipe network cell polygons from HDF file."""
        try:
            networks_path = "Geometry/Pipe Networks"
            if networks_path not in hdf_file:
                return [], []
            
            networks_group = hdf_file[networks_path]
            
            # Get network attributes
            if 'Attributes' not in networks_group:
                return [], []
            
            attributes = networks_group['Attributes'][()]
            if len(attributes) == 0:
                return [], []
            
            # Get the first network name (or could iterate through all)
            network_name = attributes[0]['Name']
            if isinstance(network_name, bytes):
                network_name = network_name.decode('utf-8', 'ignore').strip()
            
            arcpy.AddMessage(f"Processing pipe network: {network_name}")
            
            # Access the specific network group
            network_path = f"{networks_path}/{network_name}"
            if network_path not in hdf_file:
                arcpy.AddWarning(f"Network path '{network_path}' not found in HDF file")
                return [], []
            
            network_group = hdf_file[network_path]
            
            # Check for required datasets
            required_datasets = ['Cell Polygons Info', 'Cell Polygons Parts', 'Cell Polygons Points']
            for ds in required_datasets:
                if ds not in network_group:
                    arcpy.AddWarning(f"Required dataset '{ds}' not found in pipe network")
                    return [], []
            
            # Read cell polygon data
            cell_info = network_group['Cell Polygons Info'][()]
            cell_parts = network_group['Cell Polygons Parts'][()]
            cell_points = network_group['Cell Polygons Points'][()]
            
            # Read additional cell attributes if available
            cell_attributes = {}
            if 'Cell Property Table' in network_group:
                cell_property_table = network_group['Cell Property Table'][()]
                # Convert to dictionary for easier access
                for i, row in enumerate(cell_property_table):
                    cell_attributes[i] = {}
                    for field_name in cell_property_table.dtype.names:
                        value = row[field_name]
                        if isinstance(value, (bytes, np.bytes_)):
                            value = value.decode('utf-8', 'ignore').strip()
                        cell_attributes[i][field_name] = value
            
            # Read minimum elevations if available
            min_elevations = None
            if 'Cells Minimum Elevations' in network_group:
                min_elevations = network_group['Cells Minimum Elevations'][()]
            
            # Read node and conduit IDs if available
            node_conduit_ids = None
            if 'Cells Node and Conduit IDs' in network_group:
                node_conduit_ids = network_group['Cells Node and Conduit IDs'][()]
            
            valid_data, geometries = [], []
            
            # Process each cell
            for cell_idx, info in enumerate(cell_info):
                point_start_idx, point_count, part_start_idx, part_count = info
                
                try:
                    # Build polygon from parts
                    if part_count == 0:
                        continue
                    
                    parts_list = []
                    for p in range(part_start_idx, part_start_idx + part_count):
                        if p >= len(cell_parts):
                            continue
                        
                        part_info = cell_parts[p]
                        part_point_start = part_info[0]
                        part_point_count = part_info[1]
                        
                        # Extract coordinates for this part
                        coords = cell_points[part_point_start:part_point_start + part_point_count]
                        if len(coords) < 3:  # Need at least 3 points for a polygon
                            continue
                        
                        # Create arcpy array for this part
                        part_array = arcpy.Array([arcpy.Point(c[0], c[1]) for c in coords])
                        parts_list.append(part_array)
                    
                    if not parts_list:
                        continue
                    
                    # Create polygon geometry
                    if len(parts_list) == 1:
                        geom = arcpy.Polygon(parts_list[0], sr)
                    else:
                        # Multi-part polygon
                        all_parts = arcpy.Array()
                        for part in parts_list:
                            all_parts.add(part)
                        geom = arcpy.Polygon(all_parts, sr)
                    
                    # Build attribute dictionary
                    attr_dict = {
                        'cell_id': int(cell_idx),
                        'network_name': network_name
                    }
                    
                    # Add cell properties if available
                    if cell_idx in cell_attributes:
                        for key, value in cell_attributes[cell_idx].items():
                            # Clean field name
                            clean_key = key.replace(' ', '_').replace(':', '_').replace(';', '_').replace(',', '_')
                            attr_dict[clean_key] = value
                    
                    # Add minimum elevation if available
                    if min_elevations is not None and cell_idx < len(min_elevations):
                        attr_dict['min_elevation'] = float(min_elevations[cell_idx])
                    
                    # Add node and conduit IDs if available
                    if node_conduit_ids is not None and cell_idx < len(node_conduit_ids):
                        attr_dict['node_id'] = int(node_conduit_ids[cell_idx][0])
                        attr_dict['conduit_id'] = int(node_conduit_ids[cell_idx][1])
                    
                    valid_data.append(attr_dict)
                    geometries.append(geom)
                    
                except Exception as e:
                    arcpy.AddWarning(f"Error processing pipe network cell {cell_idx}: {str(e)}")
                    continue
            
            return valid_data, geometries
            
        except Exception as e:
            arcpy.AddWarning(f"Could not read pipe networks: {e}")
            return [], []

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
        
        # Get geodatabase parameters
        output_gdb = parameters[12].valueAsText
        create_gdb = parameters[13].value
        output_workspace = None
        
        # Extract project and plan info
        project_name, plan_number, base_name = extract_project_and_plan_info(hdf_path)
        
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
        
        # Setup geodatabase
        if create_gdb or output_gdb:
            if create_gdb and not output_gdb:
                # Auto-create geodatabase based on HDF name
                output_gdb = create_geodatabase_from_hdf(hdf_path, messages)
            
            # Create feature dataset with project/plan naming
            feature_dataset_name = get_feature_dataset_name(hdf_path)
            output_workspace = setup_geodatabase_output(output_gdb, feature_dataset_name, sr, messages)
            messages.addMessage(f"Output workspace set to: {output_workspace}")
        
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
                # Update output path if using geodatabase
                if output_workspace:
                    fc_name = "Breaklines2D"
                    output_fc = os.path.join(output_workspace, fc_name)
                    parameters[3].value = output_fc
                
                data, geoms = self._get_breaklines_direct(hdf_file, sr)
                fields = [("bl_id", "LONG"), ("Name", "TEXT"), ("CellSpaceNear", "FLOAT"), 
                         ("CellSpaceFar", "FLOAT"), ("NearRepeats", "LONG"), ("ProtectRadius", "LONG")]
                write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
                if output_workspace and data:
                    add_feature_class_metadata(output_fc, "2D breaklines with cell spacing attributes", hdf_path)
            
            if self.BC_LINES in geometry_elements and parameters[4].valueAsText:
                output_fc = parameters[4].valueAsText
                # Update output path if using geodatabase
                if output_workspace:
                    fc_name = "BCLines2D"
                    output_fc = os.path.join(output_workspace, fc_name)
                    parameters[4].value = output_fc
                
                data, geoms = self._get_bc_lines_direct(hdf_file, sr)
                fields = [("bc_id", "LONG"), ("Name", "TEXT"), ("Type", "TEXT")]
                write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
                if output_workspace and data:
                    add_feature_class_metadata(output_fc, "2D boundary condition lines", hdf_path)
            
            if self.PERIMETERS in geometry_elements and parameters[5].valueAsText:
                output_fc = parameters[5].valueAsText
                # Update output path if using geodatabase
                if output_workspace:
                    fc_name = "MeshPerimeters"
                    output_fc = os.path.join(output_workspace, fc_name)
                    parameters[5].value = output_fc
                
                data, geoms = self._get_mesh_areas_direct(hdf_file, sr)
                fields = [("mesh_name", "TEXT")]
                write_features_to_fc(output_fc, sr, "POLYGON", fields, data, geoms, messages)
                if output_workspace and data:
                    add_feature_class_metadata(output_fc, "2D flow area perimeter polygons", hdf_path)
            
            if self.CELL_POINTS in geometry_elements and parameters[6].valueAsText:
                output_fc = parameters[6].valueAsText
                # Update output path if using geodatabase
                if output_workspace:
                    fc_name = "MeshCellCenters"
                    output_fc = os.path.join(output_workspace, fc_name)
                    parameters[6].value = output_fc
                
                data, geoms = self._get_mesh_cell_points_direct(hdf_file, sr)
                fields = [("mesh_name", "TEXT"), ("cell_id", "LONG")]
                write_features_to_fc(output_fc, sr, "POINT", fields, data, geoms, messages)
                if output_workspace and data:
                    add_feature_class_metadata(output_fc, "Mesh cell center points", hdf_path)
            
            if self.CELL_FACES in geometry_elements and parameters[7].valueAsText:
                output_fc = parameters[7].valueAsText
                # Update output path if using geodatabase
                if output_workspace:
                    fc_name = "MeshCellFaces"
                    output_fc = os.path.join(output_workspace, fc_name)
                    parameters[7].value = output_fc
                
                if precomputed_faces:
                    data, geoms = precomputed_faces
                    fields = [("mesh_name", "TEXT"), ("face_id", "LONG")]
                    write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
                    if output_workspace and data:
                        add_feature_class_metadata(output_fc, "Mesh cell face polylines", hdf_path)
            
            if self.CELL_POLYS in geometry_elements and parameters[8].valueAsText:
                output_fc = parameters[8].valueAsText
                # Update output path if using geodatabase
                if output_workspace:
                    fc_name = "MeshCellPolygons"
                    output_fc = os.path.join(output_workspace, fc_name)
                    parameters[8].value = output_fc
                
                if precomputed_faces:
                    messages.addMessage("Constructing cell polygons from faces...")
                    data, geoms = self._get_mesh_cells_direct(hdf_file, sr, precomputed_faces, messages)
                    fields = [("mesh_name", "TEXT"), ("cell_id", "LONG")]
                    write_features_to_fc(output_fc, sr, "POLYGON", fields, data, geoms, messages)
                    if output_workspace and data:
                        add_feature_class_metadata(output_fc, "Mesh cell polygons", hdf_path)
            
            # Process pipe network elements
            # Pipe networks use the same feature dataset as other geometry
            pipe_workspace = output_workspace
            
            if self.PIPE_CONDUITS in geometry_elements and parameters[9].valueAsText:
                # Check cache first
                if not self._hdf_cache.get('has_pipe_conduits', False):
                    messages.addMessage("No pipe conduits found in the HDF file.")
                else:
                    output_fc = parameters[9].valueAsText
                    # Update output path if using geodatabase
                    if pipe_workspace:
                        fc_name = "PipeConduits"
                        output_fc = os.path.join(pipe_workspace, fc_name)
                        parameters[9].value = output_fc
                    
                    messages.addMessage("Extracting Pipe Conduits...")
                    data, geoms = self._get_pipe_conduits_direct(hdf_file, sr)
                    if data:
                        # Get dynamic fields from the data
                        fields = get_dynamic_fields_from_data(data)
                        write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
                        if pipe_workspace:
                            add_feature_class_metadata(output_fc, "Pipe conduits (storm/sewer networks)", hdf_path)
                    else:
                        messages.addMessage("No pipe conduits data extracted.")
            
            if self.PIPE_NODES in geometry_elements and parameters[10].valueAsText:
                # Check cache first
                if not self._hdf_cache.get('has_pipe_nodes', False):
                    messages.addMessage("No pipe nodes found in the HDF file.")
                else:
                    output_fc = parameters[10].valueAsText
                    # Update output path if using geodatabase
                    if pipe_workspace:
                        fc_name = "PipeNodes"
                        output_fc = os.path.join(pipe_workspace, fc_name)
                        parameters[10].value = output_fc
                    
                    messages.addMessage("Extracting Pipe Nodes...")
                    data, geoms = self._get_pipe_nodes_direct(hdf_file, sr)
                    if data:
                        # Get dynamic fields from the data
                        fields = get_dynamic_fields_from_data(data)
                        write_features_to_fc(output_fc, sr, "POINT", fields, data, geoms, messages)
                        if pipe_workspace:
                            add_feature_class_metadata(output_fc, "Pipe junction nodes", hdf_path)
                    else:
                        messages.addMessage("No pipe nodes data extracted.")
            
            if self.PIPE_NETWORKS in geometry_elements and parameters[11].valueAsText:
                # Check cache first - pipe networks might not have a specific cache flag
                output_fc = parameters[11].valueAsText
                # Update output path if using geodatabase
                if pipe_workspace:
                    fc_name = "PipeNetworks"
                    output_fc = os.path.join(pipe_workspace, fc_name)
                    parameters[11].value = output_fc
                
                messages.addMessage("Extracting Pipe Networks...")
                data, geoms = self._get_pipe_networks_direct(hdf_file, sr)
                if data:
                    # Get dynamic fields from the data
                    fields = get_dynamic_fields_from_data(data)
                    write_features_to_fc(output_fc, sr, "POLYGON", fields, data, geoms, messages)
                    if pipe_workspace:
                        add_feature_class_metadata(output_fc, "Pipe network elements", hdf_path)
                else:
                    messages.addMessage("No pipe networks data extracted.")
        
        messages.addMessage("\nProcessing complete.")
        return
    
    def getHelp(self, tool_name=None):
        """Return help documentation for the tool.
        
        This method is called when the user clicks the help button.
        It can return:
        - A URL (starting with http:// or https://)
        - A local file path (starting with file:///)
        - HTML content directly (for embedded help)
        """
        # Try local help file first
        help_file = os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 
            "Doc", "RASCommander_Help.html"
        )
        
        if os.path.exists(help_file):
            # Return local help file
            anchor = "#load-hec-ras-2d-geometry-layers"
            return f"file:///{help_file.replace(os.sep, '/')}{anchor}"
        else:
            # Fallback to online documentation
            return "https://github.com/gpt-cmdr/ras-commander-hydro#load-hec-ras-2d-geometry-layers"
    
    def getCodeSamples(self):
        """Provide code samples for using this tool programmatically."""
        return [
            {
                "title": "Basic 2D Mesh Import",
                "description": "Import mesh cell centers and perimeters",
                "code": """import arcpy

# Set input parameters
hdf_file = r"C:\\RAS_Projects\\MyProject\\MyProject.g01.hdf"
geometry_elements = ["Mesh Cell Centers", "Mesh Area Perimeters"]

# Run the tool
result = arcpy.RASCommander.LoadHECRAS2DGeometry(
    input_hdf=hdf_file,
    geometry_elements=geometry_elements,
    output_cell_points=r"memory\\MeshCenters",
    output_perimeters=r"memory\\MeshPerimeters"
)

print("2D geometry loaded successfully!")
print(f"Cell centers: {result[0]}")
print(f"Perimeters: {result[1]}")"""
            },
            {
                "title": "Complete 2D Mesh Extract",
                "description": "Extract full mesh geometry including polygons",
                "code": """import arcpy
import os

# Input HDF file
hdf_file = r"C:\\RAS_Projects\\MyProject\\MyProject.p01.hdf"

# Create geodatabase for outputs
gdb_path = os.path.join(os.path.dirname(hdf_file), "MyProject.p01.gdb")

# Extract all 2D geometry elements (including slow polygon creation)
arcpy.RASCommander.LoadHECRAS2DGeometry(
    input_hdf=hdf_file,
    geometry_elements=["2D Breaklines", "2D Boundary Condition Lines", 
                      "Mesh Area Perimeters", "Mesh Cell Centers", 
                      "Mesh Cells (Polygons)"],
    output_gdb=gdb_path,
    create_gdb=True
)

print(f"2D geometry organized in: {gdb_path}")"""
            },
            {
                "title": "Pipe Network Extraction",
                "description": "Extract storm/sewer pipe network components",
                "code": """import arcpy

# Extract pipe network elements
arcpy.RASCommander.LoadHECRAS2DGeometry(
    input_hdf=r"C:\\RAS_Projects\\Urban\\Urban.g01.hdf",
    geometry_elements=["Pipe Networks"],
    output_pipe_conduits=r"memory\\PipeConduits",
    output_pipe_nodes=r"memory\\PipeNodes"
)

# Query pipe statistics
with arcpy.da.SearchCursor("memory\\PipeConduits", ["Shape_Length", "Diameter"]) as cursor:
    total_length = 0
    for row in cursor:
        total_length += row[0]
    print(f"Total pipe length: {total_length:.2f} feet")"""
            },
            {
                "title": "Performance Optimized Extract",
                "description": "Extract mesh for large models (>100k cells)",
                "code": """import arcpy

# For large meshes, avoid polygon creation
arcpy.RASCommander.LoadHECRAS2DGeometry(
    input_hdf=r"C:\\RAS_Projects\\LargeModel.g01.hdf",
    geometry_elements=["Mesh Cell Centers", "Mesh Cell Faces", "Mesh Area Perimeters"],
    create_gdb=True
)

# Cell centers and faces provide mesh structure without expensive polygon creation
print("Large mesh geometry extracted efficiently")"""
            }
        ]
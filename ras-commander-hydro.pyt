# -*- coding: utf-8 -*-
#
# ras-commander-hydro.pyt
#
# ArcGIS Python Toolbox for HEC-RAS HDF5 Data Integration
# ===================================================================================
#
# DESCRIPTION:
# This toolbox provides tools for loading and visualizing HEC-RAS 2D geometry data
# from HDF5 files directly within ArcGIS Pro.
#
# ORIGIN AND ATTRIBUTION:
# This toolbox is a direct port of the HDF5 data extraction logic from the
# ras-commander library. All core HDF5 reading logic is derived from the
# library's HDF handling classes (e.g., HdfMesh, HdfBndry).
#
#   ras-commander library: https://github.com/gpt-cmdr/ras-commander
#
# SPONSORSHIP:
# The development and porting of these tools to the ArcGIS platform
# was generously sponsored by CLB Engineering in cooperation with ESRI.
#
#   CLB Engineering: https://clbengineering.com/
#
# ===================================================================================

import arcpy
import os
import re
import h5py
import numpy as np
from datetime import datetime, timedelta
from collections import defaultdict

# ===================================================================================
# HELPER FUNCTIONS (Dependency-free HDF logic)
# Note: These helpers use arcpy.Add* functions because they don't have access
#       to the `messages` object from the tool's `execute` method.
# ===================================================================================

def get_ras_projection_wkt(hdf_path_str: str) -> str or None:
    """
    Gets projection WKT from HDF file or an associated .prj file.
    This is a self-contained port of the HdfBase.get_projection function
    from the ras-commander library.
    """
    hdf_path = os.path.abspath(hdf_path_str)
    project_folder = os.path.dirname(hdf_path)
    wkt = None
    try:
        with h5py.File(hdf_path, 'r') as hdf_file:
            proj_wkt_attr = hdf_file.attrs.get("Projection")
            if proj_wkt_attr:
                if isinstance(proj_wkt_attr, (bytes, np.bytes_)):
                    wkt = proj_wkt_attr.decode("utf-8")
                    arcpy.AddMessage(f"Found projection in HDF file: {os.path.basename(hdf_path)}")
                    return wkt
    except Exception as e:
        arcpy.AddWarning(f"Could not read projection from HDF file attribute: {e}")
    if not wkt:
        try:
            rasmap_files = [f for f in os.listdir(project_folder) if f.lower().endswith(".rasmap")]
            if rasmap_files:
                rasmap_file_path = os.path.join(project_folder, rasmap_files[0])
                with open(rasmap_file_path, 'r', errors='ignore') as f:
                    content = f.read()
                proj_match = re.search(r'<RASProjectionFilename Filename="(.*?)"', content)
                if proj_match:
                    prj_filename = proj_match.group(1).replace('.\\', '')
                    proj_file = os.path.join(project_folder, prj_filename)
                    if os.path.exists(proj_file):
                        with open(proj_file, 'r') as f_prj:
                            wkt = f_prj.read().strip()
                            arcpy.AddMessage(f"Found projection in associated RASMapper file: {os.path.basename(proj_file)}")
                            return wkt
        except Exception as e:
            arcpy.AddWarning(f"Could not read projection from RASMapper file: {e}")
    return None


def polygonize_arcpy_optimized(line_geometries, sr):
    """
    Optimized polygon creation from line geometries using numpy arrays.
    """
    if not line_geometries:
        return None
    
    try:
        # Pre-process edges into numpy arrays for efficiency
        edge_coords = []
        edge_connections = defaultdict(list)
        tolerance = 1e-9
        
        # Extract all edge coordinates at once
        for line_idx, line in enumerate(line_geometries):
            if line is None or (hasattr(line, 'length') and line.length == 0):
                continue
            
            # Get line coordinates as numpy array
            part = line.getPart(0)
            if part.count >= 2:
                coords = np.array([[part.getObject(i).X, part.getObject(i).Y] 
                                 for i in range(part.count) if part.getObject(i)])
                
                if len(coords) >= 2:
                    edge_coords.append(coords)
                    
                    # Store connections using rounded coordinates for tolerance
                    start_key = tuple(np.round(coords[0], decimals=9))
                    end_key = tuple(np.round(coords[-1], decimals=9))
                    
                    edge_connections[start_key].append((end_key, line_idx, False))
                    edge_connections[end_key].append((start_key, line_idx, True))
        
        if not edge_coords:
            return None
        
        # Find starting point with exactly 2 connections (ideal for tracing)
        start_point = None
        for pt, connections in edge_connections.items():
            if len(connections) == 2:
                start_point = pt
                break
        
        if start_point is None:
            start_point = next(iter(edge_connections))
        
        # Trace polygon using optimized lookup
        visited = set()
        ring_coords = [start_point]
        current = start_point
        
        max_edges = len(edge_coords)
        edge_count = 0
        
        while edge_count < max_edges:
            found_next = False
            
            for next_point, edge_idx, is_reversed in edge_connections[current]:
                edge_key = (edge_idx, is_reversed)
                
                if edge_key not in visited:
                    visited.add(edge_key)
                    visited.add((edge_idx, not is_reversed))
                    
                    # Get edge coordinates
                    coords = edge_coords[edge_idx]
                    if is_reversed:
                        coords = coords[::-1]
                    
                    # Add intermediate points
                    if len(coords) > 2:
                        ring_coords.extend([tuple(c) for c in coords[1:-1]])
                    
                    ring_coords.append(next_point)
                    current = next_point
                    found_next = True
                    edge_count += 1
                    
                    # Check if closed
                    if current == start_point and len(ring_coords) > 3:
                        ring_coords = ring_coords[:-1]  # Remove duplicate
                        
                        # Convert to numpy array for efficient operations
                        ring_array = np.array(ring_coords)
                        
                        # Ensure clockwise orientation
                        if not is_clockwise_numpy(ring_array):
                            ring_array = ring_array[::-1]
                        
                        # Create polygon
                        arcpy_array = arcpy.Array([arcpy.Point(x, y) for x, y in ring_array])
                        return arcpy.Polygon(arcpy_array, sr)
                    break
            
            if not found_next:
                break
        
        # If trace failed, try to create from unique points
        if len(ring_coords) >= 3:
            # Remove duplicates while preserving order
            seen = set()
            unique_coords = []
            for coord in ring_coords:
                if coord not in seen:
                    seen.add(coord)
                    unique_coords.append(coord)
            
            if len(unique_coords) >= 3:
                ring_array = np.array(unique_coords)
                
                if not is_clockwise_numpy(ring_array):
                    ring_array = ring_array[::-1]
                
                arcpy_array = arcpy.Array([arcpy.Point(x, y) for x, y in ring_array])
                return arcpy.Polygon(arcpy_array, sr)
        
        return None
        
    except Exception as e:
        arcpy.AddWarning(f"Polygon construction failed: {e}")
        return None


def is_clockwise_numpy(coords):
    """Check if polygon coordinates are in clockwise order using numpy."""
    coords = np.asarray(coords)
    x = coords[:, 0]
    y = coords[:, 1]
    
    # Vectorized shoelace formula
    area = 0.5 * np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])
    area += 0.5 * (x[-1] * y[0] - x[0] * y[-1])
    
    return area < 0


def get_polyline_centroid_vectorized(polyline):
    """Calculate centroid of polyline using vectorized operations."""
    try:
        part = polyline.getPart(0)
        coords = np.array([[part.getObject(i).X, part.getObject(i).Y] 
                          for i in range(part.count) if part.getObject(i)])
        
        if len(coords) < 2:
            return None
        
        # Vectorized segment calculations
        segments = coords[1:] - coords[:-1]
        lengths = np.sqrt(np.sum(segments**2, axis=1))
        
        # Segment midpoints
        midpoints = (coords[:-1] + coords[1:]) / 2.0
        
        # Weighted centroid
        total_length = np.sum(lengths)
        if total_length > 0:
            weighted_coords = np.sum(midpoints * lengths[:, np.newaxis], axis=0) / total_length
            return arcpy.Point(weighted_coords[0], weighted_coords[1])
        
        return None
        
    except Exception as e:
        arcpy.AddWarning(f"Error calculating centroid: {e}")
        return None


# ===================================================================================
# ARCGIS PYTHON TOOLBOX DEFINITION
# ===================================================================================

class Toolbox(object):
    """
    ArcGIS Python Toolbox for loading HEC-RAS 2D geometry layers.
    """
    def __init__(self):
        self.label = "RAS Commander (by CLB Engineering)"
        self.alias = "RASCommander"
        self.description = "Tools for loading HEC-RAS 2D geometry from HDF5 files. Sponsored by CLB Engineering (https://clbengineering.com/)."
        self.tools = [LoadHECRAS6xHDFData]


# ===================================================================================
# TOOL: Load HEC-RAS 6.x HDF Data
# ===================================================================================

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
        self.PERIMETERS = "Mesh Area Perimeters"
        self.CELL_POINTS = "Mesh Cell Centers"
        self.CELL_FACES = "Mesh Cell Faces"
        self.CELL_POLYS = "Mesh Cells (Polygons)"
        
        # Results elements
        self.MAX_WSE_POINTS = "Max WSE at Cell Centers"
        self.MAX_FACE_VEL_POINTS = "Max Vel at Cell Faces"
        
        # Cache for HDF metadata
        self._hdf_cache = {}

    def getParameterInfo(self):
        # Separate lists for geometry and results
        geometry_elements = [self.BREAKLINES, self.PERIMETERS, self.CELL_POINTS, 
                           self.CELL_FACES, self.CELL_POLYS]
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
        params[5].value = r"memory\MeshPerimeters"
        params[6].value = r"memory\MeshCellCenters"
        params[7].value = r"memory\MeshCellFaces"
        params[8].value = r"memory\MeshCellPolygons"
        params[9].value = r"memory\MaximumWSE"
        params[10].value = r"memory\MaximumFaceVelocity"
        
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
                if field_def[1] == "DATE":
                    arcpy.management.AddField(output_fc, field_def[0], "DATE")
                else:
                    arcpy.management.AddField(output_fc, field_def[0], field_def[1], field_length=255)
                field_names.append(field_def[0])
                
        except Exception as e:
            messages.addErrorMessage(f"Failed to create output feature class {output_fc}: {e}")
            return
        
        try:
            field_names_with_shape = ["SHAPE@"] + field_names
            features_inserted = 0
            
            # Prepare all rows at once
            all_rows = []
            for i, geom in enumerate(geometries):
                if geom is None or geom.pointCount == 0:
                    continue
                
                row_data = data[i]
                row_values = [row_data.get(name) for name in field_names]
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
            
            if self.PERIMETERS in elements_to_load and parameters[5].valueAsText:
                output_fc = parameters[5].valueAsText
                data, geoms = self._get_mesh_areas_direct(hdf_file, sr)
                fields = [("mesh_name", "TEXT")]
                self._write_features_to_fc(output_fc, sr, "POLYGON", fields, data, geoms, messages)
            
            if self.CELL_POINTS in elements_to_load and parameters[6].valueAsText:
                output_fc = parameters[6].valueAsText
                data, geoms = self._get_mesh_cell_points_direct(hdf_file, sr)
                fields = [("mesh_name", "TEXT"), ("cell_id", "LONG")]
                self._write_features_to_fc(output_fc, sr, "POINT", fields, data, geoms, messages)
            
            if self.CELL_FACES in elements_to_load and parameters[7].valueAsText:
                output_fc = parameters[7].valueAsText
                if precomputed_faces:
                    data, geoms = precomputed_faces
                    fields = [("mesh_name", "TEXT"), ("face_id", "LONG")]
                    self._write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
            
            if self.CELL_POLYS in elements_to_load and parameters[8].valueAsText:
                output_fc = parameters[8].valueAsText
                if precomputed_faces:
                    messages.addMessage("Constructing cell polygons from faces...")
                    data, geoms = self._get_mesh_cells_direct(hdf_file, sr, precomputed_faces, messages)
                    fields = [("mesh_name", "TEXT"), ("cell_id", "LONG")]
                    self._write_features_to_fc(output_fc, sr, "POLYGON", fields, data, geoms, messages)
            
            # Process results elements
            if self.MAX_WSE_POINTS in elements_to_load and parameters[9].valueAsText:
                output_fc = parameters[9].valueAsText
                messages.addMessage("Extracting Maximum Water Surface Elevation points...")
                data, geoms = self._get_max_wse_points_direct(hdf_file, sr)
                fields = [("mesh_name", "TEXT"), ("cell_id", "LONG"), ("max_wse", "DOUBLE"), 
                         ("max_wse_time", "DATE")]
                self._write_features_to_fc(output_fc, sr, "POINT", fields, data, geoms, messages)
            
            if self.MAX_FACE_VEL_POINTS in elements_to_load and parameters[10].valueAsText:
                output_fc = parameters[10].valueAsText
                messages.addMessage("Extracting Maximum Face Velocity points...")
                data, geoms = self._get_max_face_velocity_points_direct(hdf_file, sr)
                fields = [("mesh_name", "TEXT"), ("face_id", "LONG"), ("max_vel", "DOUBLE"), 
                         ("time_of_max", "DATE")]
                self._write_features_to_fc(output_fc, sr, "POINT", fields, data, geoms, messages)
        
        messages.addMessage("\nProcessing complete.")
        return
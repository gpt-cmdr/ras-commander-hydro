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

import numpy as np
from collections import defaultdict

def polygonize_arcpy(line_geometries, sr):
    """
    Creates a polygon from a list of arcpy.Polyline objects using direct construction.
    This version uses numpy for efficient coordinate operations.
    """
    if not line_geometries:
        return None
    
    try:
        # Extract all edges as numpy arrays
        edges = []
        for line in line_geometries:
            if line is None or (hasattr(line, 'length') and line.length == 0):
                continue
                
            # Get the line's coordinates
            part = line.getPart(0)
            if part.count >= 2:
                # Get first and last points
                coords = []
                for i in range(part.count):
                    pt = part.getObject(i)
                    if pt:
                        coords.append([pt.X, pt.Y])
                
                if len(coords) >= 2:
                    # Store edge as numpy array
                    edge_array = np.array(coords)
                    edges.append(edge_array)
        
        if not edges:
            return None
        
        # Build connectivity graph using numpy
        # Dictionary to store point connections
        point_connections = defaultdict(list)
        tolerance = 1e-9  # Coordinate tolerance for matching points
        
        # Build adjacency information
        for edge in edges:
            if len(edge) < 2:
                continue
                
            # Get start and end points
            start = tuple(edge[0])
            end = tuple(edge[-1])
            
            # Add both directions to allow traversal
            point_connections[start].append((end, edge))
            point_connections[end].append((start, edge[::-1]))  # Reverse edge
        
        if not point_connections:
            return None
        
        # Find a starting point (preferably with exactly 2 connections)
        start_point = None
        for pt, connections in point_connections.items():
            if len(connections) == 2:
                start_point = pt
                break
        
        if start_point is None:
            # Just pick any point
            start_point = next(iter(point_connections))
        
        # Trace the polygon boundary
        ring_coords = []
        visited_edges = set()
        current = start_point
        ring_coords.append(current)
        
        max_iterations = len(edges) * 2  # Prevent infinite loops
        iterations = 0
        
        while iterations < max_iterations:
            iterations += 1
            found_next = False
            
            # Find next unvisited edge from current point
            for next_point, edge_array in point_connections[current]:
                edge_key = (current, next_point)
                if edge_key not in visited_edges:
                    visited_edges.add(edge_key)
                    visited_edges.add((next_point, current))  # Mark reverse as visited too
                    
                    # Add intermediate points from the edge (excluding start and end)
                    if len(edge_array) > 2:
                        for i in range(1, len(edge_array) - 1):
                            ring_coords.append(tuple(edge_array[i]))
                    
                    ring_coords.append(next_point)
                    current = next_point
                    found_next = True
                    
                    # Check if we've closed the ring
                    if current == start_point and len(ring_coords) > 3:
                        # Remove duplicate closing point
                        ring_coords = ring_coords[:-1]
                        
                        # Create polygon using numpy array
                        ring_array = np.array(ring_coords)
                        
                        # Ensure clockwise orientation for exterior ring
                        if not is_clockwise_numpy(ring_array):
                            ring_array = ring_array[::-1]
                        
                        # Create arcpy polygon
                        arcpy_points = [arcpy.Point(x, y) for x, y in ring_array]
                        array = arcpy.Array(arcpy_points)
                        return arcpy.Polygon(array, sr)
                    break
            
            if not found_next:
                break
        
        # If we couldn't create a closed ring, try creating from all points
        if len(ring_coords) >= 3:
            ring_array = np.array(ring_coords)
            
            # Remove duplicates while preserving order
            unique_coords = []
            seen = set()
            for coord in ring_coords:
                coord_tuple = tuple(coord)
                if coord_tuple not in seen:
                    seen.add(coord_tuple)
                    unique_coords.append(coord)
            
            if len(unique_coords) >= 3:
                ring_array = np.array(unique_coords)
                
                # Ensure clockwise orientation
                if not is_clockwise_numpy(ring_array):
                    ring_array = ring_array[::-1]
                
                arcpy_points = [arcpy.Point(x, y) for x, y in ring_array]
                array = arcpy.Array(arcpy_points)
                return arcpy.Polygon(array, sr)
        
        return None
        
    except Exception as e:
        arcpy.AddWarning(f"Direct polygon construction failed: {e}")
        return None


def is_clockwise_numpy(coords):
    """Check if polygon coordinates are in clockwise order using numpy."""
    # Use the shoelace formula to determine orientation
    # Positive area = counter-clockwise, negative = clockwise
    coords = np.array(coords)
    x = coords[:, 0]
    y = coords[:, 1]
    
    # Shoelace formula
    area = 0.5 * np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])
    area += 0.5 * (x[-1] * y[0] - x[0] * y[-1])  # Close the polygon
    
    return area < 0  # Negative area means clockwise




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
        self.tools = [
            AboutRasCommander,
            Load2DGeometryElements
        ]

# ===================================================================================
# TOOL: About RAS Commander
# ===================================================================================

class AboutRasCommander(object):
    """A tool to display information about the RAS Commander toolbox."""
    def __init__(self):
        self.label = "About RAS Commander"
        self.description = "Displays information about this toolbox, its origin, and sponsorship."
        self.canRunInBackground = False
    def getParameterInfo(self):
        return []
    def isLicensed(self):
        return True
    def execute(self, parameters, messages):
        about_message = ["===================================================================","           RAS Commander for ArcGIS","===================================================================","","This toolbox provides tools for loading HEC-RAS 2D geometry data","from HDF5 files directly within ArcGIS Pro.","","Sponsorship:","  The development of this ArcGIS toolbox was sponsored by:","  CLB Engineering (https://clbengineering.com/)","  in cooperation with ESRI.","","Origin and Attribution:","  This toolbox is based on the open-source ras-commander library.","  Find out more at: https://github.com/gpt-cmdr/ras-commander","","==================================================================="]
        for line in about_message: messages.addMessage(line)
        return

# ===================================================================================
# TOOL: Load 2D Geometry Elements
# ===================================================================================

class Load2DGeometryElements(object):
    """
    Loads one or more 2D geometry element types from a HEC-RAS HDF file.
    This tool's logic is derived from the HdfMesh and HdfBndry classes
    in the ras-commander library.
    """
    def __init__(self):
        self.label = "Load 2D Geometry Elements"
        self.description = "Extracts selected 2D geometry elements (breaklines, mesh perimeters, cells, etc.) from a HEC-RAS HDF file into separate feature classes."
        self.canRunInBackground = False
        
        self.BREAKLINES = "2D Breaklines"
        self.PERIMETERS = "Mesh Area Perimeters"
        self.CELL_POINTS = "Mesh Cell Centers"
        self.CELL_FACES = "Mesh Cell Faces"
        self.CELL_POLYS = "Mesh Cells (Polygons)"

    def getParameterInfo(self):
        element_list = [self.BREAKLINES, self.PERIMETERS, self.CELL_POINTS, self.CELL_FACES, self.CELL_POLYS]

        params = [
            arcpy.Parameter(displayName="Geometry or Plan HDF File", name="input_hdf", datatype="DEFile", parameterType="Required", direction="Input"),
            arcpy.Parameter(displayName="Geometry Elements to Load", name="elements_to_load", datatype="GPString", parameterType="Required", direction="Input", multiValue=True),
            arcpy.Parameter(displayName="Override CRS (Optional)", name="override_crs", datatype="GPSpatialReference", parameterType="Optional", direction="Input"),
            arcpy.Parameter(displayName="Output 2D Breaklines", name="output_breaklines", datatype="DEFeatureClass", parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output Mesh Area Perimeters", name="output_perimeters", datatype="DEFeatureClass", parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output Mesh Cell Centers", name="output_cell_points", datatype="DEFeatureClass", parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output Mesh Cell Faces", name="output_cell_faces", datatype="DEFeatureClass", parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output Mesh Cells (Polygons)", name="output_cell_polys", datatype="DEFeatureClass", parameterType="Optional", direction="Output", category="Outputs")
        ]
        
        params[0].filter.list = ["hdf", "g*.hdf", "p*.hdf"]
        params[1].filter.type = "ValueList"
        params[1].filter.list = element_list
        params[1].value = [self.PERIMETERS]
        
        params[3].value = r"memory\Breaklines"
        params[4].value = r"memory\MeshPerimeters"
        params[5].value = r"memory\MeshCellCenters"
        params[6].value = r"memory\MeshCellFaces"
        params[7].value = r"memory\MeshCellPolygons"
        
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    # --- HDF Data Extraction Methods ---

    def _get_breaklines_direct(self, hdf_path, sr):
            """Extracts 2D breaklines from HDF file."""
            try:
                with h5py.File(hdf_path, 'r') as hdf_file:
                    breaklines_path = "Geometry/2D Flow Area Break Lines"
                    if breaklines_path not in hdf_file: return [], []
                    bl_line_data = hdf_file[breaklines_path]
                    attributes = bl_line_data["Attributes"][()]
                    valid_data, geometries = [], []
                    for idx, (pnt_start, pnt_cnt, part_start, part_cnt) in enumerate(bl_line_data["Polyline Info"][()]):
                        attr_row = attributes[idx]
                        name = attr_row["Name"]
                        name = name.decode('utf-8', 'ignore').strip() if isinstance(name, bytes) else str(name)
                        if pnt_cnt < 2: continue
                        try:
                            points = bl_line_data["Polyline Points"][()][pnt_start:pnt_start + pnt_cnt]
                            if len(points) < 2: continue
                            if part_cnt == 1:
                                geom = arcpy.Polyline(arcpy.Array([arcpy.Point(p[0], p[1]) for p in points]), sr)
                            else:
                                parts = bl_line_data["Polyline Parts"][()][part_start:part_start + part_cnt]
                                all_parts_array = arcpy.Array()
                                for part_pnt_start, part_pnt_cnt in parts:
                                    if part_pnt_cnt > 1: all_parts_array.add(arcpy.Array([arcpy.Point(p[0], p[1]) for p in points[part_pnt_start:part_pnt_start + part_pnt_cnt]]))
                                if all_parts_array.count == 0: continue
                                geom = arcpy.Polyline(all_parts_array, sr)
                            
                            valid_data.append({
                                'bl_id': idx, 'Name': name,
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

    def _get_mesh_areas_direct(self, hdf_path, sr):
        """Extracts mesh area perimeters from HDF file."""
        try:
            with h5py.File(hdf_path, 'r') as hdf_file:
                flow_areas_path = "Geometry/2D Flow Areas"
                if flow_areas_path not in hdf_file or "Attributes" not in hdf_file[flow_areas_path]: return [], []
                attributes = hdf_file[f"{flow_areas_path}/Attributes"][()]
                mesh_area_names = [n.decode('utf-8', 'ignore').strip() for n in attributes["Name"]]
                if not mesh_area_names: return [], []
                raw_data = [{'mesh_name': name} for name in mesh_area_names]
                geometries = []
                for mesh_name in mesh_area_names:
                    perimeter_path = f"{flow_areas_path}/{mesh_name}/Perimeter"
                    if perimeter_path in hdf_file:
                        coords = hdf_file[perimeter_path][()]
                        geometries.append(arcpy.Polygon(arcpy.Array([arcpy.Point(p[0], p[1]) for p in coords]), sr))
                    else:
                        geometries.append(None)
                return raw_data, geometries
        except Exception as e:
            arcpy.AddError(f"HDF Read Error (Perimeters): {e}")
            raise arcpy.ExecuteError()
        


    def _get_mesh_cell_points_direct(self, hdf_path, sr):
        """Extracts mesh cell centers from HDF file using per-mesh approach."""
        try:
            with h5py.File(hdf_path, 'r') as hdf_file:
                flow_areas_path = "Geometry/2D Flow Areas"
                
                if flow_areas_path not in hdf_file or "Attributes" not in hdf_file[flow_areas_path]:
                    return [], []
                    
                attributes = hdf_file[f"{flow_areas_path}/Attributes"][()]
                
                # FIX: Handle the attributes array properly
                mesh_area_names = []
                num_areas = len(attributes)
                for i in range(num_areas):
                    name = attributes["Name"][i]
                    name = name.decode('utf-8', 'ignore').strip() if isinstance(name, bytes) else str(name)
                    mesh_area_names.append(name)
                
                if not mesh_area_names:
                    return [], []
                
                raw_data, geometries = [], []
                for mesh_name in mesh_area_names:
                    cell_centers_path = f"Geometry/2D Flow Areas/{mesh_name}/Cells Center Coordinate"
                    if cell_centers_path not in hdf_file:
                        arcpy.AddWarning(f"No cell center data found for mesh '{mesh_name}'")
                        continue
                        
                    cell_centers = hdf_file[cell_centers_path][()]
                    for cell_id, coords in enumerate(cell_centers):
                        raw_data.append({'mesh_name': mesh_name, 'cell_id': cell_id})
                        geometries.append(arcpy.PointGeometry(arcpy.Point(coords[0], coords[1]), sr))
                        
                return raw_data, geometries
                    
        except Exception as e:
            arcpy.AddError(f"HDF Read Error (Cell Points): {e}")
            raise arcpy.ExecuteError()


    def _get_mesh_cell_faces_direct(self, hdf_path, sr):
        """Extracts mesh cell faces from HDF file."""
        try:
            with h5py.File(hdf_path, 'r') as hdf_file:
                flow_areas_path = "Geometry/2D Flow Areas"
                if flow_areas_path not in hdf_file or "Attributes" not in hdf_file[flow_areas_path]: return [], []
                attributes = hdf_file[f"{flow_areas_path}/Attributes"][()]
                mesh_area_names = [n.decode('utf-8', 'ignore').strip() for n in attributes["Name"]]
                raw_data, geometries = [], []
                for mesh_name in mesh_area_names:
                    try:
                        base = f"Geometry/2D Flow Areas/{mesh_name}"
                        facepoints_index = hdf_file[f"{base}/Faces FacePoint Indexes"][()]
                        facepoints_coords = hdf_file[f"{base}/FacePoints Coordinate"][()]
                        faces_perim_info = hdf_file[f"{base}/Faces Perimeter Info"][()]
                        faces_perim_values = hdf_file[f"{base}/Faces Perimeter Values"][()]
                        for face_id, ((p_a, p_b), (s_row, count)) in enumerate(zip(facepoints_index, faces_perim_info)):
                            coords = [facepoints_coords[p_a]]
                            if count > 0: coords.extend(faces_perim_values[s_row : s_row + count])
                            coords.append(facepoints_coords[p_b])
                            geometries.append(arcpy.Polyline(arcpy.Array([arcpy.Point(p[0], p[1]) for p in coords]), sr))
                            raw_data.append({'mesh_name': mesh_name, 'face_id': face_id})
                    except KeyError: arcpy.AddWarning(f"No face data for mesh '{mesh_name}'.")
                return raw_data, geometries
        except Exception as e:
            arcpy.AddError(f"HDF Read Error (Cell Faces): {e}")
            raise arcpy.ExecuteError()
    


    def _get_mesh_cells_direct(self, hdf_path, sr, precomputed_faces, messages):
        """
        Optimized mesh cell extraction using numpy and direct polygon construction.
        """
        try:
            messages.addMessage("Starting optimized cell polygon creation with numpy...")
            
            with h5py.File(hdf_path, 'r') as hdf_file:
                flow_areas_path = "Geometry/2D Flow Areas"
                attributes = hdf_file[f"{flow_areas_path}/Attributes"][()]
                mesh_area_names = [n.decode('utf-8', 'ignore').strip() for n in attributes["Name"]]
                
                messages.addMessage(f"Found {len(mesh_area_names)} mesh areas: {', '.join(mesh_area_names)}")
                
                # Build face lookup with numpy arrays for coordinates
                face_lookup = {}
                face_coords_lookup = {}  # Store coordinates as numpy arrays
                
                for i, face_attr in enumerate(precomputed_faces[0]):
                    mesh_name = face_attr['mesh_name']
                    face_id = face_attr['face_id']
                    
                    if mesh_name not in face_lookup:
                        face_lookup[mesh_name] = {}
                        face_coords_lookup[mesh_name] = {}
                
                    face_geom = precomputed_faces[1][i]
                    face_lookup[mesh_name][face_id] = face_geom
                    
                    # Extract and store coordinates as numpy array
                    if face_geom and hasattr(face_geom, 'getPart'):
                        part = face_geom.getPart(0)
                        coords = []
                        for j in range(part.count):
                            pt = part.getObject(j)
                            if pt:
                                coords.append([pt.X, pt.Y])
                        if coords:
                            face_coords_lookup[mesh_name][face_id] = np.array(coords)
            
                raw_data, geometries = [], []
                total_cells_processed = 0
                
                for mesh_name in mesh_area_names:
                    messages.addMessage(f"\nProcessing mesh '{mesh_name}'...")
                    
                    try:
                        base = f"Geometry/2D Flow Areas/{mesh_name}"
                        
                        # Load cell-face relationships
                        cell_face_info = hdf_file[f"{base}/Cells Face and Orientation Info"][()]
                        cell_face_values = hdf_file[f"{base}/Cells Face and Orientation Values"][()]
                        
                        # Extract face indices and orientations as numpy arrays
                        face_indices = cell_face_values[:, 0].astype(np.int32)
                        orientations = cell_face_values[:, 1].astype(np.int32)
                        
                        mesh_faces = face_lookup.get(mesh_name, {})
                        mesh_face_coords = face_coords_lookup.get(mesh_name, {})
                        
                        cells_in_mesh = 0
                        cells_created = 0
                        
                        # Process cells in batches for progress reporting
                        batch_size = 1000
                        num_cells = len(cell_face_info)
                        
                        for batch_start in range(0, num_cells, batch_size):
                            batch_end = min(batch_start + batch_size, num_cells)
                            
                            for cell_id in range(batch_start, batch_end):
                                start, length = cell_face_info[cell_id]
                                
                                if length < 3:  # Need at least 3 faces
                                    continue
                                
                                cells_in_mesh += 1
                                
                                # Get face indices for this cell
                                cell_face_ids = face_indices[start:start + length]
                                cell_orientations = orientations[start:start + length]
                                
                                # Collect face geometries
                                face_geoms = []
                                valid_faces = True
                                
                                for j, face_id in enumerate(cell_face_ids):
                                    if face_id in mesh_faces:
                                        face_geom = mesh_faces[face_id]
                                        if face_geom:
                                            # Handle orientation if needed
                                            if cell_orientations[j] < 0 and face_id in mesh_face_coords:
                                                # Reverse the face coordinates
                                                coords = mesh_face_coords[face_id][::-1]
                                                # Create reversed polyline
                                                arcpy_points = [arcpy.Point(x, y) for x, y in coords]
                                                face_geom = arcpy.Polyline(arcpy.Array(arcpy_points), sr)
                                            face_geoms.append(face_geom)
                                    else:
                                        valid_faces = False
                                        break
                                
                                if valid_faces and len(face_geoms) >= 3:
                                    # Use optimized polygon construction
                                    polygon = polygonize_arcpy(face_geoms, sr)
                                    
                                    if polygon and polygon.area > 0:
                                        raw_data.append({'mesh_name': mesh_name, 'cell_id': cell_id})
                                        geometries.append(polygon)
                                        cells_created += 1
                                        total_cells_processed += 1
                            
                            # Progress update
                            if (batch_end % 5000) == 0 or batch_end == num_cells:
                                messages.addMessage(
                                    f"  Processed {batch_end}/{num_cells} cells in mesh '{mesh_name}' "
                                    f"({cells_created} valid polygons created)"
                                )
                        
                        messages.addMessage(
                            f"  Completed mesh '{mesh_name}': {cells_created}/{cells_in_mesh} cells created"
                        )
                        
                    except Exception as e:
                        messages.addErrorMessage(f"Error processing mesh '{mesh_name}': {str(e)}")
                        import traceback
                        messages.addErrorMessage(traceback.format_exc())
                        continue
                
                messages.addMessage(f"\nTotal cells processed: {total_cells_processed}")
                return raw_data, geometries
                
        except Exception as e:
            messages.addErrorMessage(f"Fatal error: {str(e)}")
            import traceback
            messages.addErrorMessage(traceback.format_exc())
            raise arcpy.ExecuteError()

    # --- Generic Feature Writing Method ---

    def _write_features_to_fc(self, output_fc, sr, geom_type, fields, data, geometries, messages):
        """Generic function to create and populate a feature class."""
        total_features = len(geometries)
        if total_features == 0:
            messages.addWarningMessage(f"No features found for {os.path.basename(output_fc)}. Layer will not be created.")
            return

        messages.addMessage(f"Creating feature class: {os.path.basename(output_fc)} ({total_features} features)")
        try:
            output_path, output_name = os.path.split(output_fc)
            if arcpy.Exists(output_fc): arcpy.management.Delete(output_fc)
            arcpy.management.CreateFeatureclass(output_path, output_name, geom_type, spatial_reference=sr)
            
            field_names = []
            for field_def in fields:
                arcpy.management.AddField(output_fc, field_def[0], field_def[1], field_length=255)
                field_names.append(field_def[0])
                
        except Exception as e:
            messages.addErrorMessage(f"Failed to create output feature class {output_fc}: {e}")
            return
            
        try:
            field_names_with_shape = ["SHAPE@"] + field_names
            features_inserted = 0  # Add counter
            
            with arcpy.da.InsertCursor(output_fc, field_names_with_shape) as cursor:
                for i, geom in enumerate(geometries):
                    if geom is None:
                        continue
                        
                    # Skip empty geometries
                    if geom is None or geom.pointCount == 0:
                        continue
                    
                    row_data = data[i]
                    row_values = [row_data.get(name) for name in field_names]
                    cursor.insertRow([geom] + row_values)
                    features_inserted += 1
                    
                    if features_inserted % 5000 == 0:
                        messages.addMessage(f"Inserted {features_inserted} features...")
            
            messages.addMessage(f"Successfully created {os.path.basename(output_fc)} with {features_inserted} features.")
        except Exception as e:
            messages.addErrorMessage(f"Failed during feature creation for {output_fc}: {e}")





    # --- Main Execution Logic ---
    def execute(self, parameters, messages):
        hdf_path = parameters[0].valueAsText
        elements_to_load = parameters[1].values
        
        proj_wkt = get_ras_projection_wkt(hdf_path)
        sr = None
        if proj_wkt:
            sr = arcpy.SpatialReference()
            sr.loadFromString(proj_wkt)
            messages.addMessage(f"CRS '{sr.name}' found in HEC-RAS project files.")
        elif parameters[2].value:
            sr = parameters[2].value
            messages.addMessage(f"Using user-defined override CRS: {sr.name}")
        else:
            messages.addErrorMessage("CRS could not be determined. Please use the Override CRS parameter.")
            raise arcpy.ExecuteError

        precomputed_faces = None
        if self.CELL_FACES in elements_to_load or self.CELL_POLYS in elements_to_load:
            messages.addMessage("Pre-loading Mesh Cell Faces (used by faces and polygons)...")
            precomputed_faces = self._get_mesh_cell_faces_direct(hdf_path, sr)

        if self.BREAKLINES in elements_to_load and parameters[3].valueAsText:
            output_fc = parameters[3].valueAsText
            data, geoms = self._get_breaklines_direct(hdf_path, sr)
            fields = [("bl_id", "LONG"), ("Name", "TEXT"), ("CellSpaceNear", "FLOAT"), ("CellSpaceFar", "FLOAT"), ("NearRepeats", "LONG"), ("ProtectRadius", "LONG")]
            self._write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
        
        if self.PERIMETERS in elements_to_load and parameters[4].valueAsText:
            output_fc = parameters[4].valueAsText
            data, geoms = self._get_mesh_areas_direct(hdf_path, sr)
            fields = [("mesh_name", "TEXT")]
            self._write_features_to_fc(output_fc, sr, "POLYGON", fields, data, geoms, messages)

        if self.CELL_POINTS in elements_to_load and parameters[5].valueAsText:
            output_fc = parameters[5].valueAsText
            data, geoms = self._get_mesh_cell_points_direct(hdf_path, sr)
            fields = [("mesh_name", "TEXT"), ("cell_id", "LONG")]
            self._write_features_to_fc(output_fc, sr, "POINT", fields, data, geoms, messages)

        if self.CELL_FACES in elements_to_load and parameters[6].valueAsText:
            output_fc = parameters[6].valueAsText
            if precomputed_faces:
                data, geoms = precomputed_faces
                fields = [("mesh_name", "TEXT"), ("face_id", "LONG")]
                self._write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
            else:
                messages.addWarningMessage("Could not load cell faces, pre-computation failed or was skipped.")


        if self.CELL_POLYS in elements_to_load and parameters[7].valueAsText:
            output_fc = parameters[7].valueAsText
            if precomputed_faces:
                messages.addMessage("Constructing cell polygons from faces...")
                data, geoms = self._get_mesh_cells_direct(hdf_path, sr, precomputed_faces, messages)
                fields = [("mesh_name", "TEXT"), ("cell_id", "LONG")]
                self._write_features_to_fc(output_fc, sr, "POLYGON", fields, data, geoms, messages)
            else:
                messages.addWarningMessage("Could not load cell polygons because dependent face data was not loaded.")

        messages.addMessage("\nProcessing complete.")
        return
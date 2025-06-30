# -*- coding: utf-8 -*-
"""
utils.py

Shared utility functions for RAS Commander tools.
These helpers use arcpy.Add* functions because they don't have access
to the `messages` object from the tool's `execute` method.
"""

import arcpy
import os
import re
import h5py
import numpy as np
from collections import defaultdict


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


def cache_hdf_metadata(hdf_file):
    """Pre-cache HDF metadata for faster access."""
    hdf_cache = {
        'mesh_names': [],
        'mesh_metadata': {},
        'has_results': False,
        'simulation_start_time': None
    }
    
    # Get mesh names
    flow_areas_path = "Geometry/2D Flow Areas"
    if flow_areas_path in hdf_file and "Attributes" in hdf_file[flow_areas_path]:
        attributes = hdf_file[f"{flow_areas_path}/Attributes"][()]
        hdf_cache['mesh_names'] = [n.decode('utf-8', 'ignore').strip() 
                                       for n in attributes["Name"]]
        
        # Cache mesh metadata
        for mesh_name in hdf_cache['mesh_names']:
            base_path = f"{flow_areas_path}/{mesh_name}"
            metadata = {}
            
            # Cache dataset sizes
            if f"{base_path}/Cells Center Coordinate" in hdf_file:
                metadata['cell_count'] = len(hdf_file[f"{base_path}/Cells Center Coordinate"])
            
            if f"{base_path}/Faces FacePoint Indexes" in hdf_file:
                metadata['face_count'] = len(hdf_file[f"{base_path}/Faces FacePoint Indexes"])
            
            hdf_cache['mesh_metadata'][mesh_name] = metadata
    
    # Check for boundary condition lines
    bc_lines_path = "Geometry/Boundary Condition Lines"
    if bc_lines_path in hdf_file:
        hdf_cache['has_bc_lines'] = True
        hdf_cache['bc_lines_count'] = len(hdf_file[f"{bc_lines_path}/Attributes"][()])
    else:
        hdf_cache['has_bc_lines'] = False
    
    # Check for pipe network elements
    pipe_conduits_path = "Geometry/Pipe Conduits"
    if pipe_conduits_path in hdf_file:
        hdf_cache['has_pipe_conduits'] = True
        if 'Attributes' in hdf_file[pipe_conduits_path]:
            hdf_cache['pipe_conduits_count'] = len(
                hdf_file[f"{pipe_conduits_path}/Attributes"][()]
            )
    else:
        hdf_cache['has_pipe_conduits'] = False
        
    pipe_nodes_path = "Geometry/Pipe Nodes"
    if pipe_nodes_path in hdf_file:
        hdf_cache['has_pipe_nodes'] = True
        if 'Attributes' in hdf_file[pipe_nodes_path]:
            hdf_cache['pipe_nodes_count'] = len(
                hdf_file[f"{pipe_nodes_path}/Attributes"][()]
            )
    else:
        hdf_cache['has_pipe_nodes'] = False
    
    # Check for results and get simulation time
    plan_info = hdf_file.get("Plan Data/Plan Information")
    if plan_info and 'Simulation Start Time' in plan_info.attrs:
        from datetime import datetime
        time_str = plan_info.attrs['Simulation Start Time']
        hdf_cache['simulation_start_time'] = datetime.strptime(
            time_str.decode('utf-8'), "%d%b%Y %H:%M:%S"
        )
        hdf_cache['has_results'] = True
    
    return hdf_cache


def get_dynamic_fields_from_data(data_list):
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


def write_features_to_fc(output_fc, sr, geom_type, fields, data, geometries, messages):
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
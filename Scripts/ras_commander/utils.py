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
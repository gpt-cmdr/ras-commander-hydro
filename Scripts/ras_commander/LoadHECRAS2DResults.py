# -*- coding: utf-8 -*-
"""
LoadHECRAS2DResults.py

Tool for loading HEC-RAS 2D results summary layers from HDF files including
maximum water surface elevation and face velocities.
"""

import arcpy
import os
import h5py
import numpy as np
from datetime import datetime, timedelta

# Import helper functions from utils
from .utils import (
    get_ras_projection_wkt,
    get_polyline_centroid_vectorized,
    cache_hdf_metadata,
    write_features_to_fc
)


class LoadHECRAS2DResults(object):
    """
    Loads 2D results summary data from a HEC-RAS HDF file.
    """
    def __init__(self):
        self.label = "Load HEC-RAS 2D Results Summary Layers"
        self.description = """Extracts 2D results summary data from a HEC-RAS HDF file including maximum water surface elevation and face velocities.
        
        This tool extracts summary results data from HEC-RAS plan files (p*.hdf) that contain simulation results.
        
        Available results include:
        • Max WSE at Cell Centers - Maximum water surface elevation achieved at each cell center during the simulation, with the time of occurrence
        • Max Vel at Cell Faces - Maximum velocity achieved at each cell face during the simulation, with the time of occurrence
        
        Both results types create point feature classes with attributes for the maximum value and the time when it occurred.
        
        Note: This tool requires a plan HDF file that contains results data. Geometry-only files will not work."""
        self.canRunInBackground = False
        
        # Results elements
        self.MAX_WSE_POINTS = "Max WSE at Cell Centers"
        self.MAX_FACE_VEL_POINTS = "Max Vel at Cell Faces"
        
        # Cache for HDF metadata
        self._hdf_cache = {}

    def getParameterInfo(self):
        results_elements = [self.MAX_WSE_POINTS, self.MAX_FACE_VEL_POINTS]

        params = [
            arcpy.Parameter(displayName="Plan HDF File with Results", name="input_hdf", datatype="DEFile", 
                          parameterType="Required", direction="Input"),
            arcpy.Parameter(displayName="Override CRS (Optional)", name="override_crs", datatype="GPSpatialReference", 
                          parameterType="Optional", direction="Input"),
            
            # Results elements to load
            arcpy.Parameter(displayName="Results to Load", name="results_elements", datatype="GPString", 
                          parameterType="Required", direction="Input", multiValue=True),
            
            # Output parameters
            arcpy.Parameter(displayName="Output Max WSE at Cell Centers", name="output_max_wse", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output"),
            arcpy.Parameter(displayName="Output Max Vel at Cell Faces", name="output_max_face_vel", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output")
        ]
        
        # Configure HDF file filter
        params[0].filter.list = ["hdf", "p*.hdf"]
        params[0].description = """Select a HEC-RAS plan file (p*.hdf) that contains simulation results. 
        Geometry-only files will not contain the required results data."""
        
        params[1].description = """Specify a coordinate reference system if it cannot be determined from the HEC-RAS project files. 
        The tool will first attempt to read the CRS from the HDF file or associated .prj files."""
        
        # Set filters for multi-value parameters
        params[2].filter.type = "ValueList"
        params[2].filter.list = results_elements
        params[2].value = [self.MAX_WSE_POINTS]  # Default selection
        params[2].description = """Select one or more results types to extract from the HDF file. 
        Each selected result will create a separate output feature class."""
        
        # Set default output paths and descriptions
        params[3].value = r"memory\MaximumWSE"
        params[3].description = """Output feature class for maximum water surface elevation points. 
        Includes attributes for cell ID, mesh name, maximum WSE value, and time of occurrence."""
        
        params[4].value = r"memory\MaximumFaceVelocity"
        params[4].description = """Output feature class for maximum face velocity points. 
        Includes attributes for face ID, mesh name, maximum velocity value, and time of occurrence."""
        
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
            parameters[3].enabled = self.MAX_WSE_POINTS in selected
            parameters[4].enabled = self.MAX_FACE_VEL_POINTS in selected
        return
    
    def updateMessages(self, parameters):
        """Modify messages created by internal validation."""
        # Check if HDF file is plan file with results
        if parameters[0].value:
            hdf_path = parameters[0].valueAsText
            if hdf_path and os.path.basename(hdf_path).lower().startswith('g'):
                parameters[0].setWarningMessage(
                    "This appears to be a geometry file (g*.hdf). Results data is typically in plan files (p*.hdf)."
                )
        return

    # --- HDF Data Extraction Methods ---

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

    def _get_mesh_cell_faces_direct(self, hdf_file, sr):
        """Extracts mesh cell faces needed for face velocity calculation."""
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

    # --- Main Execution Logic ---
    def execute(self, parameters, messages):
        hdf_path = parameters[0].valueAsText
        
        # Get selected elements
        results_elements = parameters[2].values if parameters[2].values else []
        
        if not results_elements:
            messages.addErrorMessage("No results elements selected for loading. Please select at least one element.")
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
            
            # Check if results exist
            if not self._hdf_cache['has_results']:
                messages.addErrorMessage("No results data found in the HDF file. Please ensure this is a plan HDF file with results.")
                raise arcpy.ExecuteError
            
            # Process results elements
            if self.MAX_WSE_POINTS in results_elements and parameters[3].valueAsText:
                output_fc = parameters[3].valueAsText
                messages.addMessage("Extracting Maximum Water Surface Elevation points...")
                data, geoms = self._get_max_wse_points_direct(hdf_file, sr)
                fields = [("mesh_name", "TEXT"), ("cell_id", "LONG"), ("max_wse", "DOUBLE"), 
                         ("max_wse_time", "DATE")]
                write_features_to_fc(output_fc, sr, "POINT", fields, data, geoms, messages)
            
            if self.MAX_FACE_VEL_POINTS in results_elements and parameters[4].valueAsText:
                output_fc = parameters[4].valueAsText
                messages.addMessage("Extracting Maximum Face Velocity points...")
                data, geoms = self._get_max_face_velocity_points_direct(hdf_file, sr)
                fields = [("mesh_name", "TEXT"), ("face_id", "LONG"), ("max_vel", "DOUBLE"), 
                         ("time_of_max", "DATE")]
                write_features_to_fc(output_fc, sr, "POINT", fields, data, geoms, messages)
        
        messages.addMessage("\nProcessing complete.")
        return
    
    def getHelp(self, tool_name):
        """Return help documentation URL for the tool."""
        help_file = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 
                                "Doc", "RASCommander_Help.html")
        if os.path.exists(help_file):
            return f"file:///{help_file.replace(os.sep, '/')}#load-hec-ras-2d-results-summary-layers"
        return None
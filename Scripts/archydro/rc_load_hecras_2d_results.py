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
from rc_utils import (
    get_ras_projection_wkt,
    get_polyline_centroid_vectorized,
    cache_hdf_metadata,
    write_features_to_fc,
    setup_geodatabase_output,
    get_unique_fc_name,
    add_feature_class_metadata,
    extract_project_and_plan_info,
    create_geodatabase_from_hdf,
    get_feature_dataset_name,
    get_feature_class_name
)


class LoadHECRAS2DResults(object):
    """
    Loads 2D results summary data from a HEC-RAS HDF file.
    """
    def __init__(self):
        # Core properties
        self.label = "Load HEC-RAS 2D Results Summary Layers"
        self.description = """Extracts 2D results summary data from a HEC-RAS HDF file including maximum water surface elevation and face velocities.
        
        This tool extracts summary results data from HEC-RAS plan files (p*.hdf) that contain simulation results.
        
        Available results include:
        • Max WSE at Cell Centers - Maximum water surface elevation achieved at each cell center during the simulation, with the time of occurrence
        • Max Vel at Cell Faces - Maximum velocity achieved at each cell face during the simulation, with the time of occurrence
        
        Both results types create point feature classes with attributes for the maximum value and the time when it occurred.
        
        Note: This tool requires a plan HDF file that contains results data. Geometry-only files will not work."""
        
        # Extended metadata properties
        self.summary = "Extract maximum WSE and velocity results from HEC-RAS 2D simulations"
        self.usage = """Select a HEC-RAS plan HDF file containing simulation results and choose which summary statistics to extract.
        
        Steps:
        1. Browse to a HEC-RAS plan file (p*.hdf) with results
        2. Select which results to extract (Max WSE, Max Velocity)
        3. Specify output locations
        4. Optionally create an organized geodatabase
        
        The tool extracts:
        • Maximum values achieved during the entire simulation
        • Time of occurrence for each maximum
        • Cell/face identification and area information
        
        Use this tool for flood mapping and hazard assessment."""
        
        # Tool behavior
        self.canRunInBackground = False
        # self.category = "HEC-RAS Results Analysis"  # REMOVE THIS LINE
        
        # Documentation and credits
        self.tags = ["HEC-RAS", "2D Results", "Water Surface Elevation", "Velocity", 
                     "Flood Mapping", "Hazard Analysis", "Arc Hydro"]
        self.credits = "CLB Engineering Corporation"
        self.author = "CLB Engineering Corporation"
        self.version = "1.0.0"
        
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
                          parameterType="Optional", direction="Output"),
            
            # Geodatabase organization parameters
            arcpy.Parameter(displayName="Output Geodatabase (Optional)", name="output_gdb", datatype="DEWorkspace", 
                          parameterType="Optional", direction="Output"),
            arcpy.Parameter(displayName="Create New Geodatabase", name="create_gdb", datatype="GPBoolean", 
                          parameterType="Optional", direction="Input")
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
        
        Attributes include:
        • Face ID and 2D area name
        • Maximum velocity magnitude (ft/s or m/s)
        • Time of maximum occurrence
        • Face location (between cell centers)
        
        Use for identifying high velocity areas and erosion potential."""
        
        # Geodatabase parameters
        params[5].description = """Specify a geodatabase to organize all output feature classes. 
        If provided, outputs will be created in this geodatabase instead of the default locations."""
        
        params[6].value = True  # Default to creating new geodatabase
        params[6].description = """Create a new geodatabase based on the HDF file name. 
        The geodatabase will be named using the pattern: ProjectName.pXX.gdb"""
        
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
        
        # Auto-populate geodatabase path when HDF file is selected
        if parameters[0].value and parameters[0].altered:  # input_hdf
            hdf_path = parameters[0].valueAsText
            
            # If create_gdb is True, auto-populate geodatabase path
            if parameters[6].value:  # create_gdb
                project_name, plan_number, base_name = extract_project_and_plan_info(hdf_path)
                gdb_name = f"{base_name}.gdb"
                gdb_path = os.path.join(os.path.dirname(hdf_path), gdb_name)
                parameters[5].value = gdb_path
        
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
        
        # Clear geodatabase validation error if create_gdb is True
        if parameters[6].value and parameters[5].hasError():  # create_gdb and output_gdb has error
            parameters[5].clearMessage()
        
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
        
        # Get geodatabase parameters
        output_gdb = parameters[5].valueAsText if len(parameters) > 5 else None
        create_gdb = parameters[6].value if len(parameters) > 6 else False
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
            
            # Check if results exist
            if not self._hdf_cache['has_results']:
                messages.addErrorMessage("No results data found in the HDF file. Please ensure this is a plan HDF file with results.")
                raise arcpy.ExecuteError
            
            # Process results elements
            if self.MAX_WSE_POINTS in results_elements and parameters[3].valueAsText:
                output_fc = parameters[3].valueAsText
                
                # Update output path if using geodatabase
                if output_workspace:
                    fc_name = "MaxWSE"
                    output_fc = os.path.join(output_workspace, fc_name)
                    parameters[3].value = output_fc
                
                messages.addMessage("Extracting Maximum Water Surface Elevation points...")
                data, geoms = self._get_max_wse_points_direct(hdf_file, sr)
                fields = [("mesh_name", "TEXT"), ("cell_id", "LONG"), ("max_wse", "DOUBLE"), 
                         ("max_wse_time", "DATE")]
                write_features_to_fc(output_fc, sr, "POINT", fields, data, geoms, messages)
                if output_workspace and data:
                    add_feature_class_metadata(output_fc, "Maximum water surface elevation at cell centers", hdf_path)
            
            if self.MAX_FACE_VEL_POINTS in results_elements and parameters[4].valueAsText:
                output_fc = parameters[4].valueAsText
                
                # Update output path if using geodatabase
                if output_workspace:
                    fc_name = "MaxVelocity"
                    output_fc = os.path.join(output_workspace, fc_name)
                    parameters[4].value = output_fc
                
                messages.addMessage("Extracting Maximum Face Velocity points...")
                data, geoms = self._get_max_face_velocity_points_direct(hdf_file, sr)
                fields = [("mesh_name", "TEXT"), ("face_id", "LONG"), ("max_vel", "DOUBLE"), 
                         ("time_of_max", "DATE")]
                write_features_to_fc(output_fc, sr, "POINT", fields, data, geoms, messages)
                if output_workspace and data:
                    add_feature_class_metadata(output_fc, "Maximum velocity at cell faces", hdf_path)
        
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
            anchor = "#load-hec-ras-2d-results-summary-layers"
            return f"file:///{help_file.replace(os.sep, '/')}{anchor}"
        else:
            # Fallback to online documentation
            return "https://github.com/gpt-cmdr/ras-commander-hydro#load-hec-ras-2d-results-summary-layers"
    
    def getCodeSamples(self):
        """Provide code samples for using this tool programmatically."""
        return [
            {
                "title": "Basic Results Extraction",
                "description": "Extract maximum WSE from simulation results",
                "code": """import arcpy

# Set input parameters
hdf_file = r"C:\\RAS_Projects\\MyProject\\MyProject.p01.hdf"
results_elements = ["Max WSE at Cell Centers"]

# Run the tool
result = arcpy.RASCommander.LoadHECRAS2DResults(
    input_hdf=hdf_file,
    results_elements=results_elements,
    output_max_wse=r"memory\\MaxWSE"
)

print("Maximum WSE extracted successfully!")

# Query statistics
with arcpy.da.SearchCursor(result[0], ["Max_WSE", "Time_of_Max"]) as cursor:
    max_wse = max(row[0] for row in cursor)
    print(f"Peak WSE in model: {max_wse:.2f} feet")"""
            },
            {
                "title": "Complete Results Analysis",
                "description": "Extract both WSE and velocity results",
                "code": """import arcpy
import os

# Input plan file with results
hdf_file = r"C:\\RAS_Projects\\MyProject\\MyProject.p01.hdf"

# Create geodatabase for results
gdb_path = os.path.join(os.path.dirname(hdf_file), "MyProject_Results.gdb")

# Extract all results
arcpy.RASCommander.LoadHECRAS2DResults(
    input_hdf=hdf_file,
    results_elements=["Max WSE at Cell Centers", "Max Vel at Cell Faces"],
    output_gdb=gdb_path,
    create_gdb=True
)

print(f"Results organized in: {gdb_path}")

# Create flood depth raster (requires Spatial Analyst)
if arcpy.CheckExtension("Spatial") == "Available":
    arcpy.CheckOutExtension("Spatial")
    wse_points = os.path.join(gdb_path, "Results", "MaxWSE_CellCenters")
    depth_raster = arcpy.sa.Idw(wse_points, "Max_WSE")
    depth_raster.save(os.path.join(gdb_path, "FloodDepth"))"""
            },
            {
                "title": "Time Analysis",
                "description": "Analyze when peak conditions occurred",
                "code": """import arcpy
from datetime import datetime

# Extract results
result = arcpy.RASCommander.LoadHECRAS2DResults(
    input_hdf=r"C:\\RAS_Projects\\TimeSeries.p01.hdf",
    results_elements=["Max WSE at Cell Centers", "Max Vel at Cell Faces"]
)

# Analyze timing of peaks
wse_fc = result[0]
vel_fc = result[1]

# Find when most cells peaked
time_counts = {}
with arcpy.da.SearchCursor(wse_fc, ["Time_of_Max"]) as cursor:
    for row in cursor:
        time_str = row[0]
        time_counts[time_str] = time_counts.get(time_str, 0) + 1

peak_time = max(time_counts, key=time_counts.get)
print(f"Most cells peaked at: {peak_time}")
print(f"Number of cells: {time_counts[peak_time]}")"""
            },
            {
                "title": "Hazard Classification",
                "description": "Classify flood hazard based on depth and velocity",
                "code": """import arcpy

# Extract results
arcpy.RASCommander.LoadHECRAS2DResults(
    input_hdf=r"C:\\RAS_Projects\\Hazard.p01.hdf",
    results_elements=["Max WSE at Cell Centers", "Max Vel at Cell Faces"],
    create_gdb=True
)

# Join velocity to WSE points for hazard analysis
# (Additional spatial join and hazard calculation code would go here)

print("Results ready for hazard classification")"""
            }
        ]
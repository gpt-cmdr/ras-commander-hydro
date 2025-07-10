# -*- coding: utf-8 -*-
"""
OrganizeRASProject.py

Master tool for organizing all HEC-RAS data from HDF files into a structured geodatabase.
This tool extracts all available geometry and results data in a single operation.
"""

import arcpy
import os
import h5py
import numpy as np

# Import helper functions from utils
from rc_utils import (
    get_ras_projection_wkt,
    setup_geodatabase_output,
    get_unique_fc_name,
    add_feature_class_metadata,
    extract_project_and_plan_info,
    create_geodatabase_from_hdf,
    get_feature_dataset_name,
    get_feature_class_name
)

# Import the individual tool classes
from rc_load_hecras_1d_geometry import LoadHECRAS1DGeometry
from rc_load_hecras_2d_geometry import LoadHECRAS2DGeometry
from rc_load_hecras_2d_results import LoadHECRAS2DResults


class OrganizeRASProject(object):
    """
    Organizes all HEC-RAS data into a structured geodatabase.
    """
    def __init__(self):
        # Core properties
        self.label = "Organize HEC-RAS Project"
        self.description = """Extracts all geometry and results from HEC-RAS files into an organized geodatabase.
        
        This tool automatically:
        • Creates a well-organized geodatabase structure
        • Extracts all available 1D geometry elements
        • Extracts all available 2D geometry elements
        • Extracts pipe network data if present
        • Extracts results data if available
        
        The geodatabase will be organized with feature datasets named by project and plan:
        • {ProjectName}_Plan{XX} - Contains all geometry and results for each plan
        
        Each plan's feature dataset will contain:
        • 1D Geometry - Cross sections, river centerlines, bank lines, structures
        • 2D Geometry - Breaklines, boundary conditions, mesh elements
        • Pipe Networks - Storm/sewer pipe networks (if present)
        • Results - Maximum WSE, velocity, and other results (if available)
        
        Note: This operation may take several minutes for large models."""
        
        # Extended metadata properties
        self.summary = "Comprehensive HEC-RAS project organization into geodatabase"
        self.usage = """Select a HEC-RAS project directory or individual HDF file to organize all data.
        
        Steps:
        1. Choose input (project folder or single HDF file)
        2. Specify output geodatabase location
        3. Select which data types to include
        4. Set processing options
        5. Tool extracts and organizes all available data
        
        Input options:
        • Project directory - Processes all HDF files found
        • Single HDF file - Processes just that file
        
        Organization structure:
        • Feature datasets by plan (Plan01, Plan02, etc.)
        • Within each plan: 1D_Geometry, 2D_Geometry, Results
        • Automatic naming and metadata
        • Preserves all HEC-RAS attributes
        
        Performance tips:
        • Disable mesh polygons for large models (>100k cells)
        • Use local drives for better performance
        • Close other applications during processing"""
        
        # Tool behavior
        self.canRunInBackground = False
        # self.category = "HEC-RAS Project Management"  # REMOVE THIS LINE
        
        # Documentation and credits
        self.tags = ["HEC-RAS", "Project Organization", "Geodatabase", "Batch Processing", 
                     "1D Geometry", "2D Geometry", "Results", "Pipe Networks", "Arc Hydro"]
        self.credits = "CLB Engineering Corporation"
        self.author = "CLB Engineering Corporation"
        self.version = "1.0.0"

    def getParameterInfo(self):
        params = [
            # Input files
            arcpy.Parameter(displayName="HEC-RAS Project Directory or HDF File", name="input_path", 
                          datatype=["DEFolder", "DEFile"], 
                          parameterType="Required", direction="Input"),
            
            # Output geodatabase
            arcpy.Parameter(displayName="Output Geodatabase", name="output_gdb", datatype="DEWorkspace", 
                          parameterType="Required", direction="Output"),
            
            # CRS override
            arcpy.Parameter(displayName="Override CRS (Optional)", name="override_crs", datatype="GPSpatialReference", 
                          parameterType="Optional", direction="Input"),
            
            # Data type selection
            arcpy.Parameter(displayName="Include 1D Geometry", name="include_1d_geometry", datatype="GPBoolean", 
                          parameterType="Optional", direction="Input"),
            arcpy.Parameter(displayName="Include 2D Geometry", name="include_2d_geometry", datatype="GPBoolean", 
                          parameterType="Optional", direction="Input"),
            arcpy.Parameter(displayName="Include 2D Results Summary Layers", name="include_2d_results", datatype="GPBoolean", 
                          parameterType="Optional", direction="Input"),
            
            # Options
            arcpy.Parameter(displayName="Include Mesh Cell Polygons", name="include_cell_polygons", datatype="GPBoolean", 
                          parameterType="Optional", direction="Input"),
        ]
        
        # Configure parameters
        params[0].description = """Select a HEC-RAS project directory to process all plan files (p*.hdf), 
        or select a single HDF file to process."""
        # params[0].category = "Input Data"  # Remove category grouping
        
        params[1].description = "Output geodatabase that will contain all extracted data in an organized structure."
        # params[1].category = "Output"  # Remove category grouping
        
        params[2].description = """Specify a coordinate reference system if it cannot be determined from the HEC-RAS project files."""
        # params[2].category = "Input Data"  # Remove category grouping
        
        # Data type selection parameters
        params[3].value = True
        params[3].description = """Include all 1D geometry elements in the output.
        
        When enabled, extracts:
        • Cross sections with station-elevation data
        • River centerlines
        • Bank lines
        • Edge lines
        • 1D structures (bridges, culverts, weirs)
        
        Disable if you only need 2D data."""
        # params[3].category = "Data Selection"  # Remove category grouping
        
        params[4].value = True
        params[4].description = """Include all 2D geometry elements in the output.
        
        When enabled, extracts:
        • 2D breaklines
        • Boundary condition lines
        • Mesh area perimeters
        • Mesh cell centers
        • Mesh cell faces
        • Mesh cell polygons (if enabled)
        • Pipe networks (conduits, nodes, and network cells)
        
        Disable if you only need 1D data."""
        # params[4].category = "Data Selection"  # Remove category grouping
        
        params[5].value = True
        params[5].description = """Include 2D results summary layers in the output.
        
        When enabled, extracts:
        • Maximum water surface elevation at cell centers
        • Maximum velocity at cell faces
        • Time of maximum occurrence
        
        Disable if you only need geometry without results."""
        # params[5].category = "Data Selection"  # Remove category grouping
        
        params[6].value = False
        params[6].description = """Include full polygon representation of mesh cells.
        
        WARNING: This can significantly increase processing time for large meshes.
        
        Performance guidelines:
        • < 10,000 cells: Safe to enable
        • 10,000 - 50,000 cells: May take several minutes
        • 50,000 - 100,000 cells: May take 10-30 minutes
        • > 100,000 cells: Not recommended (use cell centers instead)
        
        When disabled, only cell centers and faces are extracted."""
        # params[6].category = "Processing Options"  # Remove category grouping
        
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal validation."""
        # Set default output geodatabase name based on input
        if parameters[0].value and not parameters[1].altered:
            input_path = parameters[0].valueAsText
            
            if os.path.isdir(input_path):
                # For directory, use directory name
                input_name = os.path.basename(input_path.rstrip(os.sep))
                default_gdb = os.path.join(input_path, f"{input_name}_Organized.gdb")
            else:
                # For single file, use file name
                input_name = os.path.splitext(os.path.basename(input_path))[0]
                default_gdb = os.path.join(os.path.dirname(input_path), f"{input_name}_Organized.gdb")
            
            parameters[1].value = default_gdb
        
        # Enable/disable mesh polygon option based on 2D geometry selection
        if parameters[4].value == False:  # include_2d_geometry
            parameters[6].enabled = False  # include_cell_polygons
            parameters[6].value = False
        else:
            parameters[6].enabled = True
        
        return
    
    def updateMessages(self, parameters):
        """Modify messages created by internal validation."""
        # Warn about mesh polygons
        if parameters[6].value:
            parameters[6].setWarningMessage(
                "Creating cell polygons can be time-consuming for large meshes (>100,000 cells)."
            )
        
        # Warn if no data types selected
        if not any([parameters[3].value, parameters[4].value, parameters[5].value]):
            parameters[3].setErrorMessage("At least one data type must be selected.")
        
        return

    def _check_hdf_contents(self, hdf_file):
        """Check what data is available in the HDF file."""
        contents = {
            'has_1d': False,
            'has_2d': False,
            'has_pipes': False,
            'has_results': False,
            '1d_elements': [],
            '2d_elements': [],
            'pipe_elements': [],
            'result_profiles': []
        }
        
        # Check for 1D geometry
        if "Geometry/Cross Sections" in hdf_file:
            # Verify that it has actual data
            if "Geometry/Cross Sections/Attributes" in hdf_file:
                contents['has_1d'] = True
                contents['1d_elements'].append("Cross Sections")
        if "Geometry/River Centerlines" in hdf_file:
            contents['has_1d'] = True
            contents['1d_elements'].append("River Centerlines")
        if "Geometry/River Bank Lines" in hdf_file:
            contents['has_1d'] = True
            contents['1d_elements'].append("Bank Lines")
        if "Geometry/River Edge Lines" in hdf_file:
            contents['has_1d'] = True
            contents['1d_elements'].append("Edge Lines")
        if "Geometry/Structures" in hdf_file:
            # Verify that it has actual data
            if "Geometry/Structures/Attributes" in hdf_file:
                contents['has_1d'] = True
                contents['1d_elements'].append("1D Structures")
        
        # Check for 2D geometry
        if "Geometry/2D Flow Areas" in hdf_file:
            contents['has_2d'] = True
            contents['2d_elements'].append("Mesh Area Perimeters")
            contents['2d_elements'].append("Mesh Cell Centers")
            contents['2d_elements'].append("Mesh Cell Faces")
        if "Geometry/2D Flow Area Break Lines" in hdf_file:
            contents['has_2d'] = True
            contents['2d_elements'].append("2D Breaklines")
        if "Geometry/Boundary Condition Lines" in hdf_file:
            contents['has_2d'] = True
            contents['2d_elements'].append("2D Boundary Condition Lines")
        
        # Check for pipe networks - these are part of 2D geometry
        if "Geometry/Pipe Conduits" in hdf_file:
            contents['has_pipes'] = True
            contents['has_2d'] = True
            contents['pipe_elements'].append("Pipe Conduits")
            contents['2d_elements'].append("Pipe Conduits")
        if "Geometry/Pipe Nodes" in hdf_file:
            contents['has_pipes'] = True
            contents['has_2d'] = True
            contents['pipe_elements'].append("Pipe Nodes")
            contents['2d_elements'].append("Pipe Nodes")
        if "Geometry/Pipe Networks" in hdf_file:
            contents['has_pipes'] = True
            contents['has_2d'] = True
            contents['pipe_elements'].append("Pipe Networks")
            contents['2d_elements'].append("Pipe Networks")
        
        # Check for results
        if "Results/Unsteady/Output/Output Blocks/Base Output/Summary Output/2D Flow Areas" in hdf_file:
            contents['has_results'] = True
            # Get available output profiles
            results_path = "Results/Unsteady/Output/Output Blocks/Base Output/Summary Output/2D Flow Areas"
            try:
                for flow_area in hdf_file[results_path]:
                    if "Maximum Water Surface" in hdf_file[f"{results_path}/{flow_area}"]:
                        contents['result_profiles'].append("Maximum Water Surface")
                    if "Maximum Face Velocity" in hdf_file[f"{results_path}/{flow_area}"]:
                        contents['result_profiles'].append("Maximum Face Velocity")
                    break
            except:
                pass
        
        return contents

    def _load_geodatabase_to_map(self, gdb_path, messages):
        """Load all feature classes from the geodatabase into the current map, grouped by plan."""
        try:
            # Get the current project and map
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            map_obj = aprx.activeMap
            
            if not map_obj:
                messages.addWarning("No active map found. Feature classes were created but not added to the map.")
                return
            
            # Collect all feature classes and rasters first
            all_layers = []
            
            # Set workspace to the geodatabase
            arcpy.env.workspace = gdb_path
            
            # Collect all feature classes
            for dirpath, dirnames, filenames in arcpy.da.Walk(gdb_path, datatype="FeatureClass"):
                for filename in filenames:
                    fc_path = os.path.join(dirpath, filename)
                    all_layers.append((filename, fc_path, "FeatureClass"))
            
            # Collect all rasters
            for dirpath, dirnames, filenames in arcpy.da.Walk(gdb_path, datatype="RasterDataset"):
                for filename in filenames:
                    raster_path = os.path.join(dirpath, filename)
                    all_layers.append((filename, raster_path, "Raster"))
            
            # Group layers by plan
            plan_groups = {}
            for layer_name, layer_path, layer_type in all_layers:
                # Extract plan info from layer name (e.g., "ProjectName_Plan_03_LayerType")
                if "_Plan_" in layer_name:
                    parts = layer_name.split("_Plan_")
                    if len(parts) >= 2:
                        project_name = parts[0]
                        # Extract plan number (should be next 2 characters after Plan_)
                        remaining = parts[1]
                        if len(remaining) >= 2:
                            plan_number = remaining[:2]
                            plan_key = "{}_Plan_{}".format(project_name, plan_number)
                            
                            if plan_key not in plan_groups:
                                plan_groups[plan_key] = []
                            plan_groups[plan_key].append((layer_name, layer_path, layer_type))
                else:
                    # Layers without plan info go to a general group
                    if "General" not in plan_groups:
                        plan_groups["General"] = []
                    plan_groups["General"].append((layer_name, layer_path, layer_type))
            
            # Sort plans and layers within each plan
            sorted_plans = sorted(plan_groups.keys())
            
            # Add layers to map grouped by plan
            feature_classes_added = []
            messages.addMessage("\nAdding layers to map grouped by plan:")
            
            for plan_key in sorted_plans:
                # Create a group layer for this plan
                group_name = plan_key if plan_key != "General" else "Other Layers"
                
                try:
                    # Create the group layer first
                    group_layer = map_obj.createGroupLayer(group_name)
                    
                    # Sort layers within this plan alphabetically
                    plan_layers = sorted(plan_groups[plan_key], key=lambda x: x[0].lower())
                    
                    # Add layers to the group
                    layers_added_to_group = 0
                    layers_to_remove = []  # Track layers to remove after adding to group
                    
                    for layer_name, layer_path, layer_type in plan_layers:
                        try:
                            # First add the layer to the map
                            layer = map_obj.addDataFromPath(layer_path)
                            if layer:
                                # Then add it to the group layer
                                map_obj.addLayerToGroup(group_layer, layer, "AUTO_ARRANGE")
                                
                                # Track the original layer for removal
                                layers_to_remove.append(layer)
                                
                                feature_classes_added.append(layer_name)
                                layers_added_to_group += 1
                                messages.addMessage("  Added {} to {}".format(layer_name, group_name))
                        except Exception as e:
                            messages.addWarning("  Could not add {} to group: {}".format(layer_name, str(e)))
                    
                    # Remove the original layers that were added to the group
                    for layer in layers_to_remove:
                        try:
                            map_obj.removeLayer(layer)
                        except:
                            pass  # Ignore errors when removing
                    
                    if layers_added_to_group > 0:
                        messages.addMessage("Created group layer: {} with {} layers".format(group_name, layers_added_to_group))
                    else:
                        # Remove empty group layer
                        try:
                            map_obj.removeLayer(group_layer)
                        except:
                            pass
                        
                except Exception as e:
                    messages.addWarning("Could not create group layer {}: {}".format(group_name, str(e)))
                    # Fall back to adding layers to the map directly
                    plan_layers = sorted(plan_groups[plan_key], key=lambda x: x[0].lower())
                    for layer_name, layer_path, layer_type in plan_layers:
                        try:
                            layer = map_obj.addDataFromPath(layer_path)
                            if layer:
                                feature_classes_added.append(layer_name)
                                messages.addMessage("  Added {} to map (ungrouped)".format(layer_name))
                        except Exception as e:
                            messages.addWarning("  Could not add {}: {}".format(layer_name, str(e)))
            
            if feature_classes_added:
                messages.addMessage("\nSuccessfully added {} layers to the map in {} groups.".format(len(feature_classes_added), len(plan_groups)))
                
                # Try to save the project
                try:
                    aprx.save()
                    messages.addMessage("Project saved.")
                except:
                    pass  # Saving might fail in some contexts
            else:
                messages.addWarning("No feature classes were added to the map.")
                
        except Exception as e:
            messages.addWarning("Error loading geodatabase to map: {}".format(str(e)))
            messages.addWarning("Feature classes were created successfully but could not be added to the map.")

    def execute(self, parameters, messages):
        """Execute the tool."""
        input_path = parameters[0].valueAsText
        output_gdb = parameters[1].valueAsText
        override_crs = parameters[2].value
        include_1d = parameters[3].value
        include_2d = parameters[4].value
        include_results = parameters[5].value
        include_cell_polygons = parameters[6].value
        
        # Determine if input is directory or file
        hdf_files = []
        
        if os.path.isdir(input_path):
            # Find all plan files
            import glob
            pattern = os.path.join(input_path, "*.p[0-9][0-9].hdf")
            hdf_files = sorted(glob.glob(pattern))
            
            if not hdf_files:
                # Try geometry files
                pattern = os.path.join(input_path, "*.g[0-9][0-9].hdf")
                hdf_files = sorted(glob.glob(pattern))
        else:
            # Single file
            hdf_files = [input_path]
        
        if not hdf_files:
            messages.addErrorMessage("No HDF files found in the specified location.")
            return
        
        messages.addMessage(f"Found {len(hdf_files)} HDF file(s) to process")
        
        # Create output geodatabase
        gdb_folder = os.path.dirname(output_gdb)
        gdb_name = os.path.basename(output_gdb)
        
        if not arcpy.Exists(output_gdb):
            arcpy.CreateFileGDB_management(gdb_folder, gdb_name)
            messages.addMessage(f"Created geodatabase: {output_gdb}")
        
        # Process each HDF file
        for i, hdf_path in enumerate(hdf_files, 1):
            messages.addMessage(f"\n{'='*60}")
            messages.addMessage(f"Processing file {i}/{len(hdf_files)}: {os.path.basename(hdf_path)}")
            messages.addMessage(f"{'='*60}")
            
            self._process_single_hdf(hdf_path, output_gdb, override_crs, 
                                   include_1d, include_2d, include_results,
                                   include_cell_polygons, messages)
        
        messages.addMessage(f"\n{'='*60}")
        messages.addMessage(f"Processing complete! All plans organized in:")
        messages.addMessage(f"  {output_gdb}")
        
        # Load the geodatabase into the map
        messages.addMessage("\nAdding results to map...")
        self._load_geodatabase_to_map(output_gdb, messages)
    
    def _process_single_hdf(self, hdf_path, output_gdb, override_crs, 
                          include_1d, include_2d, include_results,
                          include_cell_polygons, messages):
        """Process a single HDF file."""
        # Extract project and plan info
        project_name, plan_number, base_name = extract_project_and_plan_info(hdf_path)
        
        # Get projection
        proj_wkt = get_ras_projection_wkt(hdf_path)
        sr = None
        if proj_wkt:
            sr = arcpy.SpatialReference()
            sr.loadFromString(proj_wkt)
            messages.addMessage(f"CRS '{sr.name}' found in HEC-RAS project files.")
        elif override_crs:
            sr = override_crs
            messages.addMessage(f"Using user-defined override CRS: {sr.name}")
        else:
            messages.addErrorMessage("CRS could not be determined. Please use the Override CRS parameter.")
            raise arcpy.ExecuteError
        
        # Create feature dataset for this plan
        feature_dataset_name = get_feature_dataset_name(hdf_path)
        fd_path = setup_geodatabase_output(output_gdb, feature_dataset_name, sr, messages)
        
        # Check HDF contents
        messages.addMessage("\nAnalyzing HDF file contents...")
        with h5py.File(hdf_path, 'r') as hdf_file:
            contents = self._check_hdf_contents(hdf_file)
        
        # Report what was found
        if contents['has_1d']:
            messages.addMessage(f"Found 1D geometry: {', '.join(contents['1d_elements'])}")
        if contents['has_2d']:
            messages.addMessage(f"Found 2D geometry: {', '.join(contents['2d_elements'])}")
        if contents['has_pipes']:
            messages.addMessage(f"Found pipe network data: {', '.join(contents['pipe_elements'])}")
        if contents['has_results']:
            messages.addMessage(f"Found results data with {len(set(contents['result_profiles']))} output variables")
        
        # Create mock parameters class for tool execution
        class MockParam:
            def __init__(self, value):
                self.value = value
                self.valueAsText = str(value) if value else None
                self.values = value if isinstance(value, list) else None
        
        # Process 1D Geometry
        if contents['has_1d'] and include_1d:
            messages.addMessage("\n--- Processing 1D Geometry ---")
            tool_1d = LoadHECRAS1DGeometry()
            
            # Use the plan's feature dataset
            fd_1d = fd_path
            
            # Create parameters for 1D tool
            # Note: We pass None for geodatabase parameter to prevent the individual tools
            # from overriding our carefully constructed feature class names that include
            # the project name and plan number
            params_1d = [
                hdf_path,  # input HDF
                override_crs,  # override CRS
                contents['1d_elements'],  # elements to load
                None, None, None, None, None,  # output paths (will be set individually)
                None,  # geodatabase (None to prevent tools from overriding our paths)
                False  # create_gdb (False because we already created it)
            ]
            
            mock_params = [MockParam(p) for p in params_1d]
            
            # Set output paths with full naming convention for multi-plan processing
            if "Cross Sections" in contents['1d_elements']:
                fc_name = f"{project_name}_Plan_{plan_number}_CrossSections"
                mock_params[3].valueAsText = os.path.join(fd_path, fc_name)
            if "River Centerlines" in contents['1d_elements']:
                fc_name = f"{project_name}_Plan_{plan_number}_RiverCenterlines"
                mock_params[4].valueAsText = os.path.join(fd_path, fc_name)
            if "Bank Lines" in contents['1d_elements']:
                fc_name = f"{project_name}_Plan_{plan_number}_BankLines"
                mock_params[5].valueAsText = os.path.join(fd_path, fc_name)
            if "Edge Lines" in contents['1d_elements']:
                fc_name = f"{project_name}_Plan_{plan_number}_EdgeLines"
                mock_params[6].valueAsText = os.path.join(fd_path, fc_name)
            if "1D Structures" in contents['1d_elements'] or "Hydraulic Structures" in contents['1d_elements']:
                fc_name = f"{project_name}_Plan_{plan_number}_Structures1D"
                mock_params[7].valueAsText = os.path.join(fd_path, fc_name)
            
            # Execute 1D tool
            try:
                tool_1d.execute(mock_params, messages)
            except Exception as e:
                messages.addWarning(f"Error processing 1D geometry: {e}")
                messages.addWarning("Continuing with remaining data...")
        
        # Process 2D Geometry
        if contents['has_2d'] and include_2d:
            messages.addMessage("\n--- Processing 2D Geometry ---")
            tool_2d = LoadHECRAS2DGeometry()
            
            # Determine which elements to load
            elements_2d = contents['2d_elements'][:]
            if include_cell_polygons and "Mesh Cell Centers" in elements_2d:
                elements_2d.append("Mesh Cells (Polygons)")
            
            # Create parameters for 2D tool
            params_2d = [
                hdf_path,  # input HDF
                override_crs,  # override CRS
                elements_2d,  # elements to load
                None, None, None, None, None, None, None, None, None,  # output paths
                None,  # geodatabase (None to prevent tools from overriding our paths)
                False  # create_gdb
            ]
            
            mock_params_2d = [MockParam(p) for p in params_2d]
            
            # Set output paths with full naming convention for multi-plan processing
            # Map indices to element names
            element_indices = {
                3: ("2D Breaklines", "Breaklines2D"),
                4: ("2D Boundary Condition Lines", "BCLines2D"),
                5: ("Mesh Area Perimeters", "MeshPerimeters"),
                6: ("Mesh Cell Centers", "MeshCellCenters"),
                7: ("Mesh Cell Faces", "MeshCellFaces"),
                8: ("Mesh Cells (Polygons)", "MeshCellPolygons"),
                9: ("Pipe Conduits", "PipeConduits"),
                10: ("Pipe Nodes", "PipeNodes"),
                11: ("Pipe Networks", "PipeNetworks")
            }
            
            for idx, (element_name, fc_base) in element_indices.items():
                if element_name in elements_2d:
                    fc_name = f"{project_name}_Plan_{plan_number}_{fc_base}"
                    mock_params_2d[idx].valueAsText = os.path.join(fd_path, fc_name)
            
            # Execute 2D tool
            try:
                tool_2d.execute(mock_params_2d, messages)
            except Exception as e:
                messages.addWarning(f"Error processing 2D geometry: {e}")
                messages.addWarning("Continuing with remaining data...")
        
        # Process Results
        if contents['has_results'] and include_results:
            messages.addMessage("\n--- Processing Results ---")
            tool_results = LoadHECRAS2DResults()
            
            # Create parameters for results tool
            params_results = [
                hdf_path,  # input HDF
                override_crs,  # override CRS
                ["Max WSE at Cell Centers", "Max Vel at Cell Faces"],  # results elements to load
                None,  # output max wse
                None,  # output max vel
                None,  # geodatabase (None to prevent tools from overriding our paths)
                False  # create_gdb
            ]
            
            mock_params_results = [MockParam(p) for p in params_results]
            
            # Set result output paths with full naming convention
            fc_name = f"{project_name}_Plan_{plan_number}_MaxWSE"
            mock_params_results[3].valueAsText = os.path.join(fd_path, fc_name)
            
            fc_name = f"{project_name}_Plan_{plan_number}_MaxVelocity"
            mock_params_results[4].valueAsText = os.path.join(fd_path, fc_name)
            
            # Execute results tool
            try:
                tool_results.execute(mock_params_results, messages)
            except Exception as e:
                messages.addWarning(f"Error processing results: {e}")
                messages.addWarning("Continuing with remaining data...")
        
        messages.addMessage(f"\nCompleted processing: {os.path.basename(hdf_path)}")
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
            anchor = "#organize-hec-ras-project"
            return f"file:///{help_file.replace(os.sep, '/')}{anchor}"
        else:
            # Fallback to online documentation
            return "https://github.com/gpt-cmdr/ras-commander-hydro#organize-hec-ras-project"
    
    def getCodeSamples(self):
        """Provide code samples for using this tool programmatically."""
        return [
            {
                "title": "Basic Project Organization",
                "description": "Organize a single HEC-RAS plan file",
                "code": """import arcpy
import os

# Set input HDF file
hdf_file = r"C:\\RAS_Projects\\MyProject\\MyProject.p01.hdf"

# Create output geodatabase path
gdb_path = os.path.join(os.path.dirname(hdf_file), "MyProject_Organized.gdb")

# Run the organization tool
arcpy.RASCommander.OrganizeRASProject(
    input_path=hdf_file,
    output_gdb=gdb_path,
    include_1d_geometry=True,
    include_2d_geometry=True,
    include_2d_results=True,
    include_cell_polygons=False  # Skip polygons for performance
)

print(f"Project organized in: {gdb_path}")

# List created feature datasets
arcpy.env.workspace = gdb_path
for fd in arcpy.ListDatasets("*", "Feature"):
    print(f"  Feature Dataset: {fd}")
    arcpy.env.workspace = os.path.join(gdb_path, fd)
    for fc in arcpy.ListFeatureClasses():
        print(f"    - {fc}")"""
            },
            {
                "title": "Batch Process Multiple Projects",
                "description": "Organize all HEC-RAS projects in a directory",
                "code": """import arcpy
import os

# Directory containing multiple HEC-RAS projects
projects_dir = r"C:\\RAS_Projects"

# Process each project directory
for project_name in os.listdir(projects_dir):
    project_path = os.path.join(projects_dir, project_name)
    
    if os.path.isdir(project_path):
        # Create output geodatabase
        gdb_path = os.path.join(project_path, f"{project_name}_Organized.gdb")
        
        try:
            # Organize the entire project
            arcpy.RASCommander.OrganizeRASProject(
                input_path=project_path,
                output_gdb=gdb_path,
                include_1d_geometry=True,
                include_2d_geometry=True,
                include_2d_results=True,
                include_cell_polygons=False
            )
            print(f"Organized: {project_name}")
        except Exception as e:
            print(f"Failed to organize {project_name}: {str(e)}")"""
            },
            {
                "title": "Performance Optimized",
                "description": "Organize large models efficiently",
                "code": """import arcpy

# For large models, optimize performance
arcpy.RASCommander.OrganizeRASProject(
    input_path=r"C:\\RAS_Projects\\LargeModel",
    output_gdb=r"C:\\RAS_Projects\\LargeModel_Fast.gdb",
    include_1d_geometry=True,
    include_2d_geometry=True,
    include_2d_results=True,
    include_cell_polygons=False  # Skip polygon creation for speed
)

# Post-processing tip: Create rasters from point data
# This is often faster than creating polygons
gdb = r"C:\\RAS_Projects\\LargeModel_Fast.gdb"
wse_points = os.path.join(gdb, "Plan01", "MaxWSE_CellCenters")

if arcpy.Exists(wse_points):
    # Create WSE raster (requires Spatial Analyst)
    arcpy.CheckOutExtension("Spatial")
    wse_raster = arcpy.sa.Idw(wse_points, "Max_WSE")
    wse_raster.save(os.path.join(gdb, "WSE_Raster"))"""
            },
            {
                "title": "Custom CRS Override",
                "description": "Organize with specific coordinate system",
                "code": """import arcpy

# Define custom spatial reference
custom_sr = arcpy.SpatialReference(2227)  # NAD83 State Plane CA III

# Organize with CRS override
arcpy.RASCommander.OrganizeRASProject(
    input_path=r"C:\\RAS_Projects\\StatePlane\\Project.p01.hdf",
    output_gdb=r"C:\\RAS_Projects\\StatePlane_Organized.gdb",
    override_crs=custom_sr,
    include_1d_geometry=True,
    include_2d_geometry=True,
    include_2d_results=True,
    include_cell_polygons=True
)

print("Project organized with custom coordinate system")"""
            },
            {
                "title": "Results Analysis Pipeline",
                "description": "Organize and analyze results",
                "code": """import arcpy
import os

# Organize project
project_dir = r"C:\\RAS_Projects\\FloodStudy"
gdb_path = os.path.join(project_dir, "FloodStudy_Analysis.gdb")

arcpy.RASCommander.OrganizeRASProject(
    input_path=project_dir,
    output_gdb=gdb_path,
    include_1d_geometry=True,
    include_2d_geometry=True,
    include_2d_results=True,
    include_cell_polygons=False
)

# Analyze results
arcpy.env.workspace = gdb_path
datasets = arcpy.ListDatasets("*", "Feature")

for dataset in datasets:
    wse_fc = os.path.join(gdb_path, dataset, "MaxWSE_CellCenters")
    
    if arcpy.Exists(wse_fc):
        # Get statistics
        stats = arcpy.Statistics_analysis(
            wse_fc, 
            "memory\\stats",
            [["Max_WSE", "MAX"], ["Max_WSE", "MEAN"]],
            "SA_Name"
        )
        
        print(f"\\nStatistics for {dataset}:")
        with arcpy.da.SearchCursor(stats, ["SA_Name", "MAX_Max_WSE", "MEAN_Max_WSE"]) as cursor:
            for row in cursor:
                print(f"  Area: {row[0]}")
                print(f"    Max WSE: {row[1]:.2f}")
                print(f"    Mean WSE: {row[2]:.2f}")"""
            },
            {
                "title": "2D Only Organization",
                "description": "Extract only 2D geometry and results",
                "code": """import arcpy

# Focus on 2D data only
arcpy.RASCommander.OrganizeRASProject(
    input_path=r"C:\\RAS_Projects\\2DModel",
    output_gdb=r"C:\\RAS_Projects\\2DModel_2DOnly.gdb",
    include_1d_geometry=False,  # Skip 1D
    include_2d_geometry=True,   # Include 2D geometry
    include_2d_results=True,    # Include results
    include_cell_polygons=False
)

print("2D data extracted successfully")"""
            },
            {
                "title": "Geometry Only (No Results)",
                "description": "Extract geometry without results for faster processing",
                "code": """import arcpy

# Extract geometry only - useful for model setup review
arcpy.RASCommander.OrganizeRASProject(
    input_path=r"C:\\RAS_Projects\\LargeModel.p01.hdf",
    output_gdb=r"C:\\RAS_Projects\\LargeModel_GeomOnly.gdb",
    include_1d_geometry=True,
    include_2d_geometry=True,
    include_2d_results=False,  # Skip results for speed
    include_cell_polygons=False
)

print("Geometry extracted without results")"""
            }
        ]
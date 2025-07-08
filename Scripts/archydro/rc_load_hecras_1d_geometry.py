# -*- coding: utf-8 -*-
"""
LoadHECRAS1DGeometry.py

Tool for loading HEC-RAS 1D geometry layers from HDF files including cross sections,
river centerlines, bank lines, and hydraulic structures.
"""

import arcpy
import os
import h5py
import numpy as np

# Import helper functions from utils
from rc_utils import (
    get_ras_projection_wkt,
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


class LoadHECRAS1DGeometry(object):
    """
    Loads 1D geometry elements from a HEC-RAS HDF file.
    """
    def __init__(self):
        # Core properties
        self.label = "Load HEC-RAS 1D Geometry Layers"
        self.description = """Extracts 1D geometry elements from a HEC-RAS HDF file including cross sections, river centerlines, bank lines, and hydraulic structures.
        
        This tool extracts various 1D geometry elements from HEC-RAS geometry (g*.hdf) or plan (p*.hdf) files.
        
        Available geometry elements include:
        • Cross Sections - River cross section cut lines with station-elevation data
        • River Centerlines - Main river/reach centerlines
        • Bank Lines - Left and right bank lines
        • Edge Lines - River edge lines for terrain processing
        • 1D Structures - Bridges, culverts, weirs, and other structures
        
        Note: Each selected element will create a separate feature class."""
        
        # Extended metadata properties
        self.summary = "Extract 1D geometry elements from HEC-RAS HDF files"
        self.usage = """Select a HEC-RAS geometry or plan HDF file and choose which 1D geometry elements to extract.
        
        Steps:
        1. Browse to a HEC-RAS geometry (g*.hdf) or plan (p*.hdf) file
        2. Select which geometry elements to extract
        3. Specify output locations for each selected element
        4. Optionally create an organized geodatabase
        
        The tool will automatically detect the coordinate system from the HDF file or associated .prj files."""
        
        # Tool behavior
        self.canRunInBackground = False
        self.category = "HEC-RAS Data Import"
        
        # Documentation and credits
        self.tags = ["HEC-RAS", "1D Geometry", "Cross Sections", "River Centerlines", "Hydraulic Modeling", "Arc Hydro"]
        self.credits = "CLB Engineering Corporation"
        self.author = "CLB Engineering Corporation"
        self.version = "1.0.0"
        
        # Geometry elements
        self.CROSS_SECTIONS = "Cross Sections"
        self.RIVER_CENTERLINES = "River Centerlines"
        self.BANK_LINES = "Bank Lines"
        self.EDGE_LINES = "Edge Lines"
        self.STRUCTURES = "1D Structures"
        
        # Cache for HDF metadata
        self._hdf_cache = {}

    def getParameterInfo(self):
        geometry_elements = [self.CROSS_SECTIONS, self.RIVER_CENTERLINES, self.BANK_LINES, 
                           self.EDGE_LINES, self.STRUCTURES]

        params = [
            arcpy.Parameter(displayName="Geometry or Plan HDF File", name="input_hdf", datatype="DEFile", 
                          parameterType="Required", direction="Input"),
            arcpy.Parameter(displayName="Override CRS (Optional)", name="override_crs", datatype="GPSpatialReference", 
                          parameterType="Optional", direction="Input"),
            
            # Geometry elements to load
            arcpy.Parameter(displayName="Geometry Elements to Load", name="geometry_elements", datatype="GPString", 
                          parameterType="Required", direction="Input", multiValue=True),
            
            # Output parameters
            arcpy.Parameter(displayName="Output Cross Sections", name="output_cross_sections", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output River Centerlines", name="output_centerlines", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output Bank Lines", name="output_bank_lines", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output Edge Lines", name="output_edge_lines", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs"),
            arcpy.Parameter(displayName="Output 1D Structures", name="output_structures", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs"),
            
            # Geodatabase organization parameters
            arcpy.Parameter(displayName="Output Geodatabase (Optional)", name="output_gdb", datatype="DEWorkspace", 
                          parameterType="Optional", direction="Output", category="Output"),
            arcpy.Parameter(displayName="Create New Geodatabase", name="create_gdb", datatype="GPBoolean", 
                          parameterType="Optional", direction="Input", category="Output")
        ]
        
        # Configure HDF file filter
        params[0].filter.list = ["hdf", "g*.hdf", "p*.hdf"]
        params[0].description = """Select a HEC-RAS geometry file (g*.hdf) or plan file (p*.hdf) containing 1D geometry data.
        
        The tool will automatically detect available geometry elements in the file and extract the selected ones.
        
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
        params[2].value = [self.CROSS_SECTIONS, self.RIVER_CENTERLINES]  # Default selection
        params[2].description = """Select one or more geometry elements to extract from the HDF file.
        
        Available elements:
        • Cross Sections - River cross section cut lines with detailed station-elevation data
        • River Centerlines - Main channel centerlines for each river/reach
        • Bank Lines - Left and right bank station definitions
        • Edge Lines - River edge boundaries for terrain integration
        • 1D Structures - Hydraulic structures including bridges, culverts, inline/lateral weirs
        
        Each selected element will create a separate output feature class with appropriate attributes."""
        params[2].category = "Geometry Selection"
        
        # Set default output paths and descriptions
        params[3].value = r"memory\CrossSections"
        params[3].description = """Output feature class for 1D cross sections.
        
        Attributes include:
        • River and Reach names
        • Cross section ID
        • Station locations
        • Geometry reference information"""
        
        params[4].value = r"memory\RiverCenterlines"
        params[4].description = """Output feature class for river/reach centerlines.
        
        Attributes include:
        • River name
        • Reach name
        • Length
        • Flow direction"""
        
        params[5].value = r"memory\BankLines"
        params[5].description = """Output feature class for left and right bank lines.
        
        Attributes include:
        • River and Reach names
        • Bank position (Left/Right)
        • Station references"""
        
        params[6].value = r"memory\EdgeLines"
        params[6].description = """Output feature class for river edge lines.
        
        Used for terrain modification and 2D mesh generation.
        Includes river/reach identification attributes."""
        
        params[7].value = r"memory\Structures1D"
        params[7].description = """Output feature class for 1D hydraulic structures.
        
        Structure types include:
        • Bridges
        • Culverts
        • Inline structures (weirs, gates)
        • Lateral structures
        • Pumping stations
        
        Attributes include structure type, name, and hydraulic parameters."""
        
        # Geodatabase parameters
        params[8].description = """Specify a geodatabase to organize all output feature classes.
        
        If provided:
        • All outputs will be created in feature datasets within this geodatabase
        • Feature datasets will be organized by geometry type
        • Automatic naming conventions will be applied
        
        Leave empty to use default output locations."""
        params[8].category = "Output Organization"
        
        params[9].value = True  # Default to creating new geodatabase
        params[9].description = """Create a new geodatabase based on the HDF file name.
        
        When enabled:
        • Creates geodatabase named: ProjectName.pXX.gdb
        • Organizes outputs in feature datasets
        • Maintains HEC-RAS project structure
        • Preserves all attribute relationships
        
        Recommended for organizing complete HEC-RAS projects."""
        params[9].category = "Output Organization"
        
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
            parameters[3].enabled = self.CROSS_SECTIONS in selected
            parameters[4].enabled = self.RIVER_CENTERLINES in selected
            parameters[5].enabled = self.BANK_LINES in selected
            parameters[6].enabled = self.EDGE_LINES in selected
            parameters[7].enabled = self.STRUCTURES in selected
            
        # Auto-populate geodatabase path when HDF file is selected
        if parameters[0].value and parameters[0].altered:  # input_hdf
            hdf_path = parameters[0].valueAsText
            
            # If create_gdb is True, auto-populate geodatabase path
            if parameters[9].value:  # create_gdb
                project_name, plan_number, base_name = extract_project_and_plan_info(hdf_path)
                gdb_name = f"{base_name}.gdb"
                gdb_path = os.path.join(os.path.dirname(hdf_path), gdb_name)
                parameters[8].value = gdb_path
        return
    
    def updateMessages(self, parameters):
        """Modify messages created by internal validation."""
        # Clear geodatabase validation error if create_gdb is True
        if parameters[9].value and parameters[8].hasError():  # create_gdb and output_gdb has error
            parameters[8].clearMessage()
        return

    # --- HDF Data Extraction Methods ---

    def _get_cross_sections_direct(self, hdf_file, sr):
        """Extracts cross sections from HDF file."""
        try:
            xs_path = "Geometry/Cross Sections"
            if xs_path not in hdf_file:
                arcpy.AddMessage("No cross sections found in HDF file.")
                return [], []
            
            # Check if required datasets exist
            required_datasets = ["Attributes", "Polyline Info", "Polyline Points", 
                               "Station Elevation Info", "Station Elevation Values"]
            for dataset in required_datasets:
                if f"{xs_path}/{dataset}" not in hdf_file:
                    arcpy.AddWarning(f"Cross sections data incomplete: missing '{dataset}' dataset.")
                    return [], []
            
            # Get attributes
            attributes = hdf_file[f"{xs_path}/Attributes"][()]
            
            # Get polyline geometry
            polyline_info = hdf_file[f"{xs_path}/Polyline Info"][()]
            polyline_points = hdf_file[f"{xs_path}/Polyline Points"][()]
            
            # Get station-elevation data
            sta_elev_info = hdf_file[f"{xs_path}/Station Elevation Info"][()]
            sta_elev_values = hdf_file[f"{xs_path}/Station Elevation Values"][()]
            
            # Get Manning's n data
            mannings_info = hdf_file[f"{xs_path}/Manning's n Info"][()]
            mannings_values = hdf_file[f"{xs_path}/Manning's n Values"][()]
            
            valid_data, geometries = [], []
            
            for idx, attr in enumerate(attributes):
                # Get polyline info
                pnt_start, pnt_cnt, _, _ = polyline_info[idx]
                
                if pnt_cnt < 2:
                    continue
                
                # Extract points and create polyline
                points = polyline_points[pnt_start:pnt_start + pnt_cnt]
                arcpy_array = arcpy.Array([arcpy.Point(p[0], p[1]) for p in points])
                geom = arcpy.Polyline(arcpy_array, sr)
                
                # Extract attributes
                river = attr["River"].decode('utf-8', 'ignore').strip()
                reach = attr["Reach"].decode('utf-8', 'ignore').strip()
                rs = attr["RS"].decode('utf-8', 'ignore').strip()
                
                # Get station-elevation profile
                se_start, se_count = sta_elev_info[idx]
                sta_elev_pairs = []
                if se_count > 0:
                    se_data = sta_elev_values[se_start:se_start + se_count]
                    sta_elev_pairs = [(float(s), float(e)) for s, e in se_data]
                
                # Get Manning's n profile
                mn_start, mn_count = mannings_info[idx]
                mannings_pairs = []
                if mn_count > 0:
                    mn_data = mannings_values[mn_start:mn_start + mn_count]
                    mannings_pairs = [(float(s), float(n)) for s, n in mn_data]
                
                valid_data.append({
                    'xs_id': int(idx),
                    'River': river,
                    'Reach': reach,
                    'RS': rs,
                    'LeftBank': float(attr["Left Bank"]),
                    'RightBank': float(attr["Right Bank"]),
                    'LenLeft': float(attr["Len Left"]),
                    'LenChannel': float(attr["Len Channel"]),
                    'LenRight': float(attr["Len Right"]),
                    'StationElevation': str(sta_elev_pairs)[:255],  # Convert to string for field storage
                    'ManningsN': str(mannings_pairs)[:255]
                })
                geometries.append(geom)
            
            return valid_data, geometries
            
        except Exception as e:
            arcpy.AddWarning(f"Error reading cross sections: {e}")
            return [], []

    def _get_river_centerlines_direct(self, hdf_file, sr):
        """Extracts river centerlines from HDF file."""
        try:
            centerlines_path = "Geometry/River Centerlines"
            if centerlines_path not in hdf_file:
                return [], []
            
            # Get attributes
            attributes = hdf_file[f"{centerlines_path}/Attributes"][()]
            
            # Get polyline geometry
            polyline_info = hdf_file[f"{centerlines_path}/Polyline Info"][()]
            polyline_points = hdf_file[f"{centerlines_path}/Polyline Points"][()]
            
            valid_data, geometries = [], []
            
            for idx, attr in enumerate(attributes):
                # Get polyline info
                pnt_start, pnt_cnt, _, _ = polyline_info[idx]
                
                if pnt_cnt < 2:
                    continue
                
                # Extract points and create polyline
                points = polyline_points[pnt_start:pnt_start + pnt_cnt]
                arcpy_array = arcpy.Array([arcpy.Point(p[0], p[1]) for p in points])
                geom = arcpy.Polyline(arcpy_array, sr)
                
                # Extract attributes
                river_name = attr["River Name"].decode('utf-8', 'ignore').strip()
                reach_name = attr["Reach Name"].decode('utf-8', 'ignore').strip()
                
                valid_data.append({
                    'river_id': int(idx),
                    'RiverName': river_name,
                    'ReachName': reach_name,
                    'USType': attr["US Type"].decode('utf-8', 'ignore').strip(),
                    'DSType': attr["DS Type"].decode('utf-8', 'ignore').strip()
                })
                geometries.append(geom)
            
            return valid_data, geometries
            
        except Exception as e:
            arcpy.AddError(f"HDF Read Error (River Centerlines): {e}")
            raise arcpy.ExecuteError("Failed to read river centerlines from HDF file")

    def _get_bank_lines_direct(self, hdf_file, sr):
        """Extracts bank lines from HDF file."""
        try:
            bank_lines_path = "Geometry/River Bank Lines"
            if bank_lines_path not in hdf_file:
                arcpy.AddMessage("No bank lines found in HDF file.")
                return [], []
            
            # Get polyline geometry
            polyline_info = hdf_file[f"{bank_lines_path}/Polyline Info"][()]
            polyline_points = hdf_file[f"{bank_lines_path}/Polyline Points"][()]
            
            valid_data, geometries = [], []
            
            # Bank lines typically come in pairs (left and right)
            bank_sides = ['Left', 'Right']
            
            for idx in range(len(polyline_info)):
                # Get polyline info
                pnt_start, pnt_cnt, _, _ = polyline_info[idx]
                
                if pnt_cnt < 2:
                    continue
                
                # Extract points and create polyline
                points = polyline_points[pnt_start:pnt_start + pnt_cnt]
                arcpy_array = arcpy.Array([arcpy.Point(p[0], p[1]) for p in points])
                geom = arcpy.Polyline(arcpy_array, sr)
                
                # Determine bank side
                bank_side = bank_sides[idx % 2] if idx < len(bank_sides) else f"Bank_{idx}"
                
                valid_data.append({
                    'bank_id': int(idx),
                    'BankSide': bank_side,
                    'Length': float(geom.length)
                })
                geometries.append(geom)
            
            return valid_data, geometries
            
        except Exception as e:
            arcpy.AddWarning(f"Error reading bank lines: {e}")
            return [], []

    def _get_edge_lines_direct(self, hdf_file, sr):
        """Extracts edge lines from HDF file."""
        try:
            edge_lines_path = "Geometry/River Edge Lines"
            if edge_lines_path not in hdf_file:
                arcpy.AddMessage("No edge lines found in HDF file.")
                return [], []
            
            # Get polyline geometry
            polyline_info = hdf_file[f"{edge_lines_path}/Polyline Info"][()]
            polyline_points = hdf_file[f"{edge_lines_path}/Polyline Points"][()]
            
            valid_data, geometries = [], []
            
            for idx in range(len(polyline_info)):
                # Get polyline info
                pnt_start, pnt_cnt, _, _ = polyline_info[idx]
                
                if pnt_cnt < 2:
                    continue
                
                # Extract points and create polyline
                points = polyline_points[pnt_start:pnt_start + pnt_cnt]
                arcpy_array = arcpy.Array([arcpy.Point(p[0], p[1]) for p in points])
                geom = arcpy.Polyline(arcpy_array, sr)
                
                valid_data.append({
                    'edge_id': int(idx),
                    'EdgeType': f"Edge_{idx}",
                    'Length': float(geom.length)
                })
                geometries.append(geom)
            
            return valid_data, geometries
            
        except Exception as e:
            arcpy.AddWarning(f"Error reading edge lines: {e}")
            return [], []

    def _get_structures_direct(self, hdf_file, sr):
        """Extracts hydraulic structures from HDF file."""
        try:
            structures_path = "Geometry/Structures"
            if structures_path not in hdf_file:
                arcpy.AddMessage("No hydraulic structures found in HDF file.")
                return [], []
            
            # Get attributes
            attributes = hdf_file[f"{structures_path}/Attributes"][()]
            
            # Get centerline geometry
            centerline_info = hdf_file[f"{structures_path}/Centerline Info"][()]
            centerline_points = hdf_file[f"{structures_path}/Centerline Points"][()]
            
            valid_data, geometries = [], []
            
            for idx, attr in enumerate(attributes):
                # Get centerline info
                pnt_start, pnt_cnt, _, _ = centerline_info[idx]
                
                if pnt_cnt < 2:
                    continue
                
                # Extract points and create polyline
                points = centerline_points[pnt_start:pnt_start + pnt_cnt]
                arcpy_array = arcpy.Array([arcpy.Point(p[0], p[1]) for p in points])
                geom = arcpy.Polyline(arcpy_array, sr)
                
                # Extract attributes
                struct_type = attr["Type"].decode('utf-8', 'ignore').strip()
                river = attr["River"].decode('utf-8', 'ignore').strip()
                reach = attr["Reach"].decode('utf-8', 'ignore').strip()
                rs = attr["RS"].decode('utf-8', 'ignore').strip()
                
                valid_data.append({
                    'struct_id': int(idx),
                    'Type': struct_type,
                    'River': river,
                    'Reach': reach,
                    'RS': rs,
                    'Description': attr["Description"].decode('utf-8', 'ignore').strip()[:255]
                })
                geometries.append(geom)
            
            return valid_data, geometries
            
        except Exception as e:
            arcpy.AddWarning(f"Error reading 1D structures: {e}")
            return [], []

    # --- Main Execution Logic ---
    def execute(self, parameters, messages):
        hdf_path = parameters[0].valueAsText
        
        # Get selected elements
        geometry_elements = parameters[2].values if parameters[2].values else []
        
        if not geometry_elements:
            messages.addErrorMessage("No geometry elements selected for loading. Please select at least one element.")
            raise arcpy.ExecuteError
        
        # Get geodatabase parameters
        output_gdb = parameters[8].valueAsText
        create_gdb = parameters[9].value
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
        
        # Open HDF file once
        with h5py.File(hdf_path, 'r') as hdf_file:
            messages.addMessage("Reading HDF file structure...")
            
            # Process geometry elements
            if self.CROSS_SECTIONS in geometry_elements and parameters[3].valueAsText:
                output_fc = parameters[3].valueAsText
                # Update output path if using geodatabase
                if output_workspace:
                    fc_name = "CrossSections"
                    output_fc = os.path.join(output_workspace, fc_name)
                    parameters[3].value = output_fc
                
                messages.addMessage("Extracting Cross Sections...")
                data, geoms = self._get_cross_sections_direct(hdf_file, sr)
                fields = [("xs_id", "LONG"), ("River", "TEXT"), ("Reach", "TEXT"), 
                         ("RS", "TEXT"), ("LeftBank", "DOUBLE"), ("RightBank", "DOUBLE"),
                         ("LenLeft", "DOUBLE"), ("LenChannel", "DOUBLE"), ("LenRight", "DOUBLE"),
                         ("StationElevation", "TEXT", 255), ("ManningsN", "TEXT", 255)]
                write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
                if output_workspace and data:
                    add_feature_class_metadata(output_fc, "River cross section cut lines with station-elevation data", hdf_path)
            
            if self.RIVER_CENTERLINES in geometry_elements and parameters[4].valueAsText:
                output_fc = parameters[4].valueAsText
                # Update output path if using geodatabase
                if output_workspace:
                    fc_name = "RiverCenterlines"
                    output_fc = os.path.join(output_workspace, fc_name)
                    parameters[4].value = output_fc
                
                messages.addMessage("Extracting River Centerlines...")
                data, geoms = self._get_river_centerlines_direct(hdf_file, sr)
                fields = [("river_id", "LONG"), ("RiverName", "TEXT"), ("ReachName", "TEXT"),
                         ("USType", "TEXT"), ("DSType", "TEXT")]
                write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
                if output_workspace and data:
                    add_feature_class_metadata(output_fc, "Main river/reach centerlines", hdf_path)
            
            if self.BANK_LINES in geometry_elements and parameters[5].valueAsText:
                output_fc = parameters[5].valueAsText
                # Update output path if using geodatabase
                if output_workspace:
                    fc_name = "BankLines"
                    output_fc = os.path.join(output_workspace, fc_name)
                    parameters[5].value = output_fc
                
                messages.addMessage("Extracting Bank Lines...")
                data, geoms = self._get_bank_lines_direct(hdf_file, sr)
                fields = [("bank_id", "LONG"), ("BankSide", "TEXT"), ("Length", "DOUBLE")]
                write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
                if output_workspace and data:
                    add_feature_class_metadata(output_fc, "Left and right bank lines", hdf_path)
            
            if self.EDGE_LINES in geometry_elements and parameters[6].valueAsText:
                output_fc = parameters[6].valueAsText
                # Update output path if using geodatabase
                if output_workspace:
                    fc_name = "EdgeLines"
                    output_fc = os.path.join(output_workspace, fc_name)
                    parameters[6].value = output_fc
                
                messages.addMessage("Extracting Edge Lines...")
                data, geoms = self._get_edge_lines_direct(hdf_file, sr)
                fields = [("edge_id", "LONG"), ("EdgeType", "TEXT"), ("Length", "DOUBLE")]
                write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
                if output_workspace and data:
                    add_feature_class_metadata(output_fc, "River edge lines for terrain processing", hdf_path)
            
            if self.STRUCTURES in geometry_elements and parameters[7].valueAsText:
                output_fc = parameters[7].valueAsText
                # Update output path if using geodatabase
                if output_workspace:
                    fc_name = "Structures1D"
                    output_fc = os.path.join(output_workspace, fc_name)
                    parameters[7].value = output_fc
                
                messages.addMessage("Extracting 1D Structures...")
                data, geoms = self._get_structures_direct(hdf_file, sr)
                fields = [("struct_id", "LONG"), ("Type", "TEXT"), ("River", "TEXT"),
                         ("Reach", "TEXT"), ("RS", "TEXT"), ("Description", "TEXT", 255)]
                write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
                if output_workspace and data:
                    add_feature_class_metadata(output_fc, "Bridges, culverts, weirs, and other 1D structures", hdf_path)
        
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
            anchor = "#load-hec-ras-1d-geometry-layers"
            return f"file:///{help_file.replace(os.sep, '/')}{anchor}"
        else:
            # Fallback to online documentation
            return "https://github.com/gpt-cmdr/ras-commander-hydro#load-hec-ras-1d-geometry-layers"
    
    def getCodeSamples(self):
        """Provide code samples for using this tool programmatically."""
        return [
            {
                "title": "Basic 1D Geometry Import",
                "description": "Import cross sections and river centerlines to memory",
                "code": """import arcpy

# Set input parameters
hdf_file = r"C:\\RAS_Projects\\MyProject\\MyProject.g01.hdf"
geometry_elements = ["Cross Sections", "River Centerlines"]

# Run the tool
result = arcpy.RASCommander.LoadHECRAS1DGeometry(
    input_hdf=hdf_file,
    geometry_elements=geometry_elements,
    output_cross_sections=r"memory\\CrossSections",
    output_centerlines=r"memory\\RiverCenterlines"
)

print("1D geometry loaded successfully!")
print(f"Cross sections: {result[0]}")
print(f"Centerlines: {result[1]}")"""
            },
            {
                "title": "Organize to Geodatabase",
                "description": "Extract all 1D geometry to an organized geodatabase",
                "code": """import arcpy
import os

# Input HDF file
hdf_file = r"C:\\RAS_Projects\\MyProject\\MyProject.p01.hdf"

# Create geodatabase path
gdb_path = os.path.join(os.path.dirname(hdf_file), "MyProject.p01.gdb")

# Extract all 1D geometry elements
arcpy.RASCommander.LoadHECRAS1DGeometry(
    input_hdf=hdf_file,
    geometry_elements=["Cross Sections", "River Centerlines", "Bank Lines", "1D Structures"],
    output_gdb=gdb_path,
    create_gdb=True
)

print(f"1D geometry organized in: {gdb_path}")"""
            },
            {
                "title": "With Custom Projection",
                "description": "Load geometry with a specific coordinate system",
                "code": """import arcpy

# Define custom spatial reference
sr = arcpy.SpatialReference(26915)  # NAD83 UTM Zone 15N

# Run tool with override CRS
arcpy.RASCommander.LoadHECRAS1DGeometry(
    input_hdf=r"C:\\RAS_Projects\\MyProject.g01.hdf",
    override_crs=sr,
    geometry_elements=["Cross Sections", "River Centerlines"],
    create_gdb=True
)"""
            }
        ]
Looking at the existing toolbox structure and the QGIS plugin examples, I'll create a new tool for extracting 1D HEC-RAS geometry data. Here's the implementation:

## 1. First, create the new tool file: `Scripts/ras_commander/LoadHECRAS1DGeometry.py`

```python
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
from .utils import (
    get_ras_projection_wkt,
    cache_hdf_metadata,
    write_features_to_fc,
    get_dynamic_fields_from_data
)


class LoadHECRAS1DGeometry(object):
    """
    Loads 1D geometry elements from a HEC-RAS HDF file.
    """
    def __init__(self):
        self.label = "Load HEC-RAS 1D Geometry Layers"
        self.description = """Extracts 1D geometry elements from a HEC-RAS HDF file including cross sections, river centerlines, bank lines, and hydraulic structures.
        
        This tool extracts various 1D geometry elements from HEC-RAS geometry (g*.hdf) or plan (p*.hdf) files.
        
        Available geometry elements include:
        • Cross Sections - River cross section cut lines with station-elevation data
        • River Centerlines - Main river/reach centerlines
        • Bank Lines - Left and right bank lines
        • Edge Lines - River edge lines for terrain processing
        • Hydraulic Structures - Bridges, culverts, weirs, and other structures
        
        Note: Each selected element will create a separate feature class."""
        self.canRunInBackground = False
        
        # Geometry elements
        self.CROSS_SECTIONS = "Cross Sections"
        self.RIVER_CENTERLINES = "River Centerlines"
        self.BANK_LINES = "Bank Lines"
        self.EDGE_LINES = "Edge Lines"
        self.STRUCTURES = "Hydraulic Structures"
        
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
            arcpy.Parameter(displayName="Output Hydraulic Structures", name="output_structures", datatype="DEFeatureClass", 
                          parameterType="Optional", direction="Output", category="Outputs")
        ]
        
        # Configure HDF file filter
        params[0].filter.list = ["hdf", "g*.hdf", "p*.hdf"]
        params[0].description = "Select a HEC-RAS geometry file (g*.hdf) or plan file (p*.hdf) containing 1D geometry data."
        
        params[1].description = """Specify a coordinate reference system if it cannot be determined from the HEC-RAS project files. 
        The tool will first attempt to read the CRS from the HDF file or associated .prj files."""
        
        # Set filters for multi-value parameters
        params[2].filter.type = "ValueList"
        params[2].filter.list = geometry_elements
        params[2].value = [self.CROSS_SECTIONS, self.RIVER_CENTERLINES]  # Default selection
        params[2].description = """Select one or more geometry elements to extract from the HDF file. 
        Each selected element will create a separate output feature class."""
        
        # Set default output paths and descriptions
        params[3].value = r"memory\CrossSections"
        params[3].description = "Output feature class for 1D cross sections with attributes."
        
        params[4].value = r"memory\RiverCenterlines"
        params[4].description = "Output feature class for river/reach centerlines."
        
        params[5].value = r"memory\BankLines"
        params[5].description = "Output feature class for left and right bank lines."
        
        params[6].value = r"memory\EdgeLines"
        params[6].description = "Output feature class for river edge lines."
        
        params[7].value = r"memory\HydraulicStructures"
        params[7].description = "Output feature class for hydraulic structures (bridges, culverts, etc.)."
        
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
        return
    
    def updateMessages(self, parameters):
        """Modify messages created by internal validation."""
        return

    # --- HDF Data Extraction Methods ---

    def _get_cross_sections_direct(self, hdf_file, sr):
        """Extracts cross sections from HDF file."""
        try:
            xs_path = "Geometry/Cross Sections"
            if xs_path not in hdf_file:
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
            arcpy.AddError(f"HDF Read Error (Cross Sections): {e}")
            raise arcpy.ExecuteError("Failed to read cross sections from HDF file")

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
            arcpy.AddError(f"HDF Read Error (Bank Lines): {e}")
            raise arcpy.ExecuteError("Failed to read bank lines from HDF file")

    def _get_edge_lines_direct(self, hdf_file, sr):
        """Extracts edge lines from HDF file."""
        try:
            edge_lines_path = "Geometry/River Edge Lines"
            if edge_lines_path not in hdf_file:
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
            arcpy.AddError(f"HDF Read Error (Edge Lines): {e}")
            raise arcpy.ExecuteError("Failed to read edge lines from HDF file")

    def _get_structures_direct(self, hdf_file, sr):
        """Extracts hydraulic structures from HDF file."""
        try:
            structures_path = "Geometry/Structures"
            if structures_path not in hdf_file:
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
            arcpy.AddError(f"HDF Read Error (Structures): {e}")
            raise arcpy.ExecuteError("Failed to read structures from HDF file")

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
        
        # Open HDF file once
        with h5py.File(hdf_path, 'r') as hdf_file:
            messages.addMessage("Reading HDF file structure...")
            
            # Process geometry elements
            if self.CROSS_SECTIONS in geometry_elements and parameters[3].valueAsText:
                output_fc = parameters[3].valueAsText
                messages.addMessage("Extracting Cross Sections...")
                data, geoms = self._get_cross_sections_direct(hdf_file, sr)
                fields = [("xs_id", "LONG"), ("River", "TEXT"), ("Reach", "TEXT"), 
                         ("RS", "TEXT"), ("LeftBank", "DOUBLE"), ("RightBank", "DOUBLE"),
                         ("LenLeft", "DOUBLE"), ("LenChannel", "DOUBLE"), ("LenRight", "DOUBLE"),
                         ("StationElevation", "TEXT", 255), ("ManningsN", "TEXT", 255)]
                write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
            
            if self.RIVER_CENTERLINES in geometry_elements and parameters[4].valueAsText:
                output_fc = parameters[4].valueAsText
                messages.addMessage("Extracting River Centerlines...")
                data, geoms = self._get_river_centerlines_direct(hdf_file, sr)
                fields = [("river_id", "LONG"), ("RiverName", "TEXT"), ("ReachName", "TEXT"),
                         ("USType", "TEXT"), ("DSType", "TEXT")]
                write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
            
            if self.BANK_LINES in geometry_elements and parameters[5].valueAsText:
                output_fc = parameters[5].valueAsText
                messages.addMessage("Extracting Bank Lines...")
                data, geoms = self._get_bank_lines_direct(hdf_file, sr)
                fields = [("bank_id", "LONG"), ("BankSide", "TEXT"), ("Length", "DOUBLE")]
                write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
            
            if self.EDGE_LINES in geometry_elements and parameters[6].valueAsText:
                output_fc = parameters[6].valueAsText
                messages.addMessage("Extracting Edge Lines...")
                data, geoms = self._get_edge_lines_direct(hdf_file, sr)
                fields = [("edge_id", "LONG"), ("EdgeType", "TEXT"), ("Length", "DOUBLE")]
                write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
            
            if self.STRUCTURES in geometry_elements and parameters[7].valueAsText:
                output_fc = parameters[7].valueAsText
                messages.addMessage("Extracting Hydraulic Structures...")
                data, geoms = self._get_structures_direct(hdf_file, sr)
                fields = [("struct_id", "LONG"), ("Type", "TEXT"), ("River", "TEXT"),
                         ("Reach", "TEXT"), ("RS", "TEXT"), ("Description", "TEXT", 255)]
                write_features_to_fc(output_fc, sr, "POLYLINE", fields, data, geoms, messages)
        
        messages.addMessage("\nProcessing complete.")
        return
    
    def getHelp(self, tool_name):
        """Return help documentation URL for the tool."""
        help_file = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 
                                "Doc", "RASCommander_Help.html")
        if os.path.exists(help_file):
            return f"file:///{help_file.replace(os.sep, '/')}#load-hec-ras-1d-geometry-layers"
        return None
```

## 2. Update `Scripts/ras_commander/__init__.py`

```python
# -*- coding: utf-8 -*-
"""
RAS Commander tools package.

This package contains the tools for loading HEC-RAS data into ArcGIS Pro.
"""

from .LoadRASTerrain import LoadRASTerrain
from .LoadHECRAS2DGeometry import LoadHECRAS2DGeometry
from .LoadHECRAS2DResults import LoadHECRAS2DResults
from .LoadHECRAS1DGeometry import LoadHECRAS1DGeometry

__all__ = ['LoadRASTerrain', 'LoadHECRAS2DGeometry', 'LoadHECRAS2DResults', 'LoadHECRAS1DGeometry']
```

## 3. Update `toolboxes/RAS Commander.pyt`

```python
# -*- coding: utf-8 -*-
#
# RAS Commander.pyt
#
# ArcGIS Python Toolbox for HEC-RAS HDF5 Data Integration
# ===================================================================================
#
# DESCRIPTION:
# This toolbox provides tools for loading and visualizing HEC-RAS 1D and 2D geometry,
# terrain, and results data from HDF5 files directly within ArcGIS Pro.
#
# ORIGIN AND ATTRIBUTION:
# This toolbox is a direct port of the HDF5 data extraction logic from the
# ras-commander library. All core HDF5 reading logic is derived from the
# library's HDF handling classes (e.g., HdfMesh, HdfBndry, HdfXsec, HdfStruc).
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

import sys
import os

# Add the Scripts directory to the Python path so we can import our modules
toolbox_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.join(os.path.dirname(toolbox_dir), 'Scripts')
if scripts_dir not in sys.path:
    sys.path.insert(0, scripts_dir)

# Import the tool classes from our modules
from ras_commander import LoadRASTerrain, LoadHECRAS2DGeometry, LoadHECRAS2DResults, LoadHECRAS1DGeometry


class Toolbox(object):
    """
    ArcGIS Python Toolbox for loading HEC-RAS 1D and 2D geometry, terrain, and results layers.
    """
    def __init__(self):
        self.label = "RAS Commander"
        self.alias = "RASCommander"
        self.description = "Tools for loading HEC-RAS 1D and 2D geometry, terrain, and results from HDF5 files. Sponsored by CLB Engineering (https://clbengineering.com/)."
        # List the tool classes
        self.tools = [LoadHECRAS1DGeometry, LoadHECRAS2DGeometry, LoadHECRAS2DResults, LoadRASTerrain]
```

## 4. Update the Help Documentation

Add a new section to `Doc/RASCommander_Help.html` after the Overview section:

```html
<div class="tool-section">
    <h2 id="load-hec-ras-1d-geometry-layers">Load HEC-RAS 1D Geometry Layers</h2>
    <p>This tool extracts 1D geometry elements from HEC-RAS geometry (<code>g*.hdf</code>) or plan (<code>p*.hdf</code>) files.</p>
    
    <h3>Available Geometry Elements</h3>
    <ul>
        <li><strong>Cross Sections</strong> - River cross section cut lines with station-elevation profiles</li>
        <li><strong>River Centerlines</strong> - Main river and reach centerlines</li>
        <li><strong>Bank Lines</strong> - Left and right bank lines</li>
        <li><strong>Edge Lines</strong> - River edge lines for terrain processing</li>
        <li><strong>Hydraulic Structures</strong> - Bridges, culverts, weirs, and other structures</li>
    </ul>
    
    <div class="note">
        <h4>Performance Note</h4>
        <p>Cross sections with many station-elevation points may take longer to process. The station-elevation and Manning's n profiles are stored as text fields due to their variable length.</p>
    </div>
    
    <h3>Parameters</h3>
    <div class="parameter">
        <h4>Geometry or Plan HDF File</h4>
        <p>Select a HEC-RAS geometry file (<code>g*.hdf</code>) or plan file (<code>p*.hdf</code>) containing 1D geometry data.</p>
    </div>
    <div class="parameter">
        <h4>Override CRS (Optional)</h4>
        <p>Specify a coordinate reference system if it cannot be determined from the HEC-RAS project files. The tool will first attempt to read the CRS from the HDF file or associated .prj files.</p>
    </div>
    <div class="parameter">
        <h4>Geometry Elements to Load</h4>
        <p>Select one or more geometry elements to extract. Each selected element will create a separate output feature class.</p>
    </div>
</div>
```

## 5. Create XML metadata file: `toolboxes/RAS Commander.LoadHECRAS1DGeometry.pyt.xml`

```xml
<?xml version="1.0" encoding="UTF-8"?>
<metadata xml:lang="en">
    <Esri>
        <CreaDate>20250629</CreaDate>
        <CreaTime>14000000</CreaTime>
        <ArcGISFormat>1.0</ArcGISFormat>
        <SyncOnce>TRUE</SyncOnce>
        <ModDate>20250629</ModDate>
        <ModTime>140000</ModTime>
    </Esri>
    <tool name="LoadHECRAS1DGeometry" displayname="Load HEC-RAS 1D Geometry Layers" toolboxalias="RASCommander">
        <summary>Extracts 1D geometry elements from a HEC-RAS HDF file including cross sections, river centerlines, bank lines, and hydraulic structures.</summary>
        <arcToolboxHelpPath>toolboxes\RAS Commander.pyt</arcToolboxHelpPath>
        <parameters>
            <param name="input_hdf" displayname="Geometry or Plan HDF File" datatype="File" direction="Input" expression="input_hdf" type="Required">
                <dialogReference>Select a HEC-RAS geometry file (g*.hdf) or plan file (p*.hdf) containing 1D geometry data.</dialogReference>
            </param>
            <param name="override_crs" displayname="Override CRS (Optional)" datatype="Spatial Reference" direction="Input" expression="override_crs" type="Optional">
                <dialogReference>Specify a coordinate reference system if it cannot be determined from the HEC-RAS project files. The tool will first attempt to read the CRS from the HDF file or associated .prj files.</dialogReference>
            </param>
            <param name="geometry_elements" displayname="Geometry Elements to Load" datatype="String" direction="Input" expression="geometry_elements" type="Required">
                <dialogReference>Select one or more geometry elements to extract from the HDF file. Each selected element will create a separate output feature class.</dialogReference>
            </param>
        </parameters>
    </tool>
    <dataIdInfo>
        <idCitation>
            <resTitle>Load HEC-RAS 1D Geometry Layers</resTitle>
        </idCitation>
        <idAbs>This tool extracts various 1D geometry elements from HEC-RAS geometry (g*.hdf) or plan (p*.hdf) files.

Available geometry elements include:
• Cross Sections - River cross section cut lines with station-elevation data
• River Centerlines - Main river/reach centerlines  
• Bank Lines - Left and right bank lines
• Edge Lines - River edge lines for terrain processing
• Hydraulic Structures - Bridges, culverts, weirs, and other structures

Select one or more elements to extract. Each selected element will create a separate feature class.

Note: Station-elevation and Manning's n profiles are stored as text fields.</idAbs>
        <idPurp>To extract and visualize HEC-RAS 1D geometry elements in ArcGIS Pro for analysis and mapping.</idPurp>
        <searchKeys>
            <keyword>HEC-RAS</keyword>
            <keyword>1D Geometry</keyword>
            <keyword>Cross Sections</keyword>
            <keyword>River Centerlines</keyword>
            <keyword>Bank Lines</keyword>
            <keyword>Hydraulic Structures</keyword>
        </searchKeys>
    </dataIdInfo>
</metadata>
```

This implementation provides a comprehensive tool for extracting 1D HEC-RAS geometry data that:

1. **Follows the existing toolbox patterns** - Similar structure to LoadHECRAS2DGeometry
2. **Reads directly from HDF files** - No dependency on external ras-commander library functions
3. **Handles multiple geometry types** - Cross sections, centerlines, bank lines, edge lines, and structures
4. **Preserves important attributes** - Including station-elevation profiles and Manning's n values
5. **Provides flexible output options** - Users can select which elements to extract
6. **Includes proper error handling** - Graceful handling of missing data
7. **Supports CRS detection and override** - Same pattern as other tools

The tool integrates seamlessly with the existing RAS Commander toolbox and provides access to all major 1D geometry elements from HEC-RAS models.
# -*- coding: utf-8 -*-
"""
LoadRASTerrain.py

Tool for loading HEC-RAS terrain layers from a project's rasmap file.
"""

import arcpy
import os
import xml.etree.ElementTree as ET
from pathlib import Path


class LoadRASTerrain(object):
    """
    Loads one or more HEC-RAS terrain layers from a project's rasmap file.
    
    IMPORTANT: This tool loads the underlying terrain TIFFs as a VRT (Virtual Raster) 
    with the priority from HEC-RAS in place. It does NOT include terrain modifications 
    done as vector terrain modifications in RAS Mapper. Only the base terrain raster 
    data will be loaded.
    """
    def __init__(self):
        # Core properties
        self.label = "Load HEC-RAS Terrain"
        self.description = """Loads terrain layers defined in a HEC-RAS project's .rasmap file into the current map.
        
        ⚠️ IMPORTANT REMINDER: The loaded terrain layers are base VRT files only. Any vector terrain 
        modifications (breaklines, high ground, etc.) made in RAS Mapper are NOT included in these layers."""
        
        # Extended metadata properties
        self.summary = "Import HEC-RAS terrain layers from RAS Mapper VRT files"
        self.usage = """Select a HEC-RAS project file to load associated terrain layers.
        
        Steps:
        1. Browse to a HEC-RAS project file (*.prj)
        2. Choose to import all terrains or select specific ones
        3. Terrain layers will be added to the current map
        
        Limitations:
        • Only loads base terrain rasters (VRT files)
        • Does NOT include vector terrain modifications
        • Does NOT include breaklines or high ground modifications
        • For complete terrain with modifications, use HEC-RAS itself
        
        The tool reads terrain definitions from the .rasmap file and loads
        the corresponding Virtual Raster (VRT) files."""
        
        # Tool behavior
        self.canRunInBackground = False
        self.category = "HEC-RAS Data Import"
        
        # Documentation and credits
        self.tags = ["HEC-RAS", "Terrain", "DEM", "RAS Mapper", "VRT", "Elevation", "Arc Hydro"]
        self.credits = "CLB Engineering Corporation"
        self.author = "CLB Engineering Corporation"
        self.version = "1.0.0"
        
        # Cache for terrain metadata
        self._terrain_cache = {}

    def getParameterInfo(self):
        params = [
            arcpy.Parameter(
                displayName="HEC-RAS Project File (*.prj)",
                name="in_ras_project",
                datatype="DEFile",
                parameterType="Required",
                direction="Input"
            ),
            arcpy.Parameter(
                displayName="Import All Terrains",
                name="import_all",
                datatype="GPBoolean",
                parameterType="Optional",
                direction="Input"
            ),
            arcpy.Parameter(
                displayName="Terrains to Load",
                name="terrains_to_load",
                datatype="GPString",
                parameterType="Optional",
                direction="Input",
                multiValue=True
            )
        ]
        
        params[0].filter.list = ["prj"]
        params[0].description = """Select the HEC-RAS project file (*.prj).
        
        The tool will:
        1. Automatically find the associated .rasmap file
        2. Read terrain layer definitions
        3. Locate the corresponding VRT files
        
        Ensure the project has been opened in RAS Mapper at least once
        to generate the necessary terrain files."""
        params[0].category = "Input Data"
        
        params[1].value = False
        params[1].description = """Check this box to import all terrain layers found in the project.
        
        When enabled:
        • All terrain layers will be loaded
        • Terrain selection list will be disabled
        • Useful for complete project visualization
        
        When disabled:
        • Select specific terrain layers to load
        • Reduces memory usage for large projects"""
        params[1].category = "Import Options"
        
        params[2].filter.type = "ValueList"
        params[2].description = """Select specific terrain layers to load.
        
        Available terrains:
        • Listed from the .rasmap file
        • Each represents a VRT (Virtual Raster)
        • May include primary terrain and alternates
        
        Note: Disabled when 'Import All Terrains' is checked.
        
        ⚠️ Remember: These are base terrain files only.
        Vector modifications are not included."""
        params[2].category = "Import Options"
        
        return params
    
    def updateMessages(self, parameters):
        """Modify the messages created by internal parameter validation."""
        return

    def isLicensed(self):
        return True

    def _get_rasmap_path_from_prj(self, prj_path_str):
        """Finds the .rasmap file associated with a .prj file."""
        if not prj_path_str or not os.path.exists(prj_path_str):
            return None
        
        p = Path(prj_path_str)
        rasmap_path = p.with_suffix('.rasmap')
        
        if rasmap_path.exists():
            return str(rasmap_path)
        
        arcpy.AddWarning(f"Could not find associated .rasmap file for {p.name}")
        return None

    def _get_terrain_info_from_rasmap(self, rasmap_path_str):
        """Parses a .rasmap file to get terrain names and their HDF file paths."""
        if not rasmap_path_str or not os.path.exists(rasmap_path_str):
            return {}
            
        terrains = {}
        try:
            tree = ET.parse(rasmap_path_str)
            root = tree.getroot()
            project_folder = os.path.dirname(rasmap_path_str)
            
            for layer in root.findall(".//Terrains/Layer"):
                name = layer.get('Name')
                filename_rel = layer.get('Filename')
                
                if name and filename_rel:
                    # HEC-RAS paths can start with '.\', remove it for robust joining
                    clean_rel_path = filename_rel.lstrip('.\\/')
                    filename_abs = os.path.join(project_folder, clean_rel_path)
                    terrains[name] = os.path.normpath(filename_abs)
                    
        except ET.ParseError as e:
            arcpy.AddWarning(f"Error parsing .rasmap file {os.path.basename(rasmap_path_str)}: {e}")
        except Exception as e:
            arcpy.AddWarning(f"An unexpected error occurred while reading the .rasmap file: {e}")
            
        return terrains

    def updateParameters(self, parameters):
        """Modify the parameters on the GUI according to user input."""
        if parameters[0].value and parameters[0].altered:
            prj_path = parameters[0].valueAsText
            rasmap_path = self._get_rasmap_path_from_prj(prj_path)
            
            if rasmap_path:
                self._terrain_cache = self._get_terrain_info_from_rasmap(rasmap_path)
                terrain_names = list(self._terrain_cache.keys())
                parameters[2].filter.list = sorted(terrain_names)
            else:
                parameters[2].filter.list = []
                self._terrain_cache = {}

        if parameters[1].value is True:
            parameters[2].enabled = False
            parameters[2].value = None # Clear selection if "Import All" is checked
        else:
            parameters[2].enabled = True
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        prj_path = parameters[0].valueAsText
        import_all = parameters[1].value
        selected_terrains = parameters[2].values

        if not self._terrain_cache:
            rasmap_path = self._get_rasmap_path_from_prj(prj_path)
            if rasmap_path:
                self._terrain_cache = self._get_terrain_info_from_rasmap(rasmap_path)

        if not self._terrain_cache:
            messages.addErrorMessage("No terrains found or .rasmap file could not be read. Aborting.")
            return

        terrains_to_load = []
        if import_all:
            terrains_to_load = list(self._terrain_cache.keys())
            messages.addMessage("Import All selected. Loading all available terrains...")
        elif selected_terrains:
            terrains_to_load = selected_terrains
            messages.addMessage(f"Loading selected terrains: {', '.join(terrains_to_load)}")
        else:
            messages.addErrorMessage("No terrains selected for loading.")
            return
            
        try:
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            # Ensure there is an active map
            if not aprx.activeMap:
                 messages.addErrorMessage("No active map found. Please open a map view and try again.")
                 return
            map = aprx.activeMap
        except Exception as e:
            messages.addErrorMessage(f"Could not access the current ArcGIS Pro project or map: {e}")
            return
            
        layers_added = 0
        for terrain_name in terrains_to_load:
            hdf_path_str = self._terrain_cache.get(terrain_name)
            if not hdf_path_str:
                messages.addWarningMessage(f"Could not find path for terrain '{terrain_name}'. Skipping.")
                continue
            
            p = Path(hdf_path_str)
            vrt_path = str(p.with_suffix('.vrt'))
            
            if os.path.exists(vrt_path):
                try:
                    map.addDataFromPath(vrt_path)
                    messages.addMessage(f"Successfully added terrain layer: {terrain_name}")
                    layers_added += 1
                except Exception as e:
                    messages.addWarningMessage(f"Failed to add layer for '{terrain_name}' from path {vrt_path}: {e}")
            else:
                messages.addWarningMessage(f"Associated VRT file not found for terrain '{terrain_name}'. Expected at: {vrt_path}")
        
        messages.addMessage(f"\nProcessing complete. Added {layers_added} terrain layer(s).")
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
            anchor = "#load-hec-ras-terrain"
            return f"file:///{help_file.replace(os.sep, '/')}{anchor}"
        else:
            # Fallback to online documentation
            return "https://github.com/gpt-cmdr/ras-commander-hydro#load-hec-ras-terrain"
    
    def getCodeSamples(self):
        """Provide code samples for using this tool programmatically."""
        return [
            {
                "title": "Load All Terrains",
                "description": "Import all terrain layers from a HEC-RAS project",
                "code": """import arcpy

# Set input project file
project_file = r"C:\\RAS_Projects\\MyProject\\MyProject.prj"

# Load all terrains
arcpy.RASCommander.LoadRASTerrain(
    in_ras_project=project_file,
    import_all=True
)

print("All terrain layers loaded to current map")

# List loaded raster layers
aprx = arcpy.mp.ArcGISProject("CURRENT")
map_obj = aprx.activeMap
for lyr in map_obj.listLayers():
    if lyr.isRasterLayer:
        print(f"  - {lyr.name}")"""
            },
            {
                "title": "Load Specific Terrains",
                "description": "Select and load specific terrain layers",
                "code": """import arcpy

# Set input project file
project_file = r"C:\\RAS_Projects\\MyProject\\MyProject.prj"

# Load specific terrains
terrains = ["Primary Terrain", "Proposed Conditions"]

arcpy.RASCommander.LoadRASTerrain(
    in_ras_project=project_file,
    import_all=False,
    terrains_to_load=terrains
)

print(f"Loaded {len(terrains)} terrain layers")"""
            },
            {
                "title": "Terrain Analysis",
                "description": "Load terrain and perform elevation analysis",
                "code": """import arcpy

# Load primary terrain
project_file = r"C:\\RAS_Projects\\FloodStudy\\FloodStudy.prj"

arcpy.RASCommander.LoadRASTerrain(
    in_ras_project=project_file,
    import_all=False,
    terrains_to_load=["Primary Terrain"]
)

# Get loaded terrain layer
aprx = arcpy.mp.ArcGISProject("CURRENT")
map_obj = aprx.activeMap
terrain_lyr = None

for lyr in map_obj.listLayers():
    if lyr.isRasterLayer and "Primary Terrain" in lyr.name:
        terrain_lyr = lyr
        break

if terrain_lyr:
    # Get elevation statistics
    desc = arcpy.Describe(terrain_lyr)
    print(f"Terrain extent: {desc.extent}")
    
    # Calculate statistics if needed
    arcpy.management.CalculateStatistics(terrain_lyr)"""
            },
            {
                "title": "Batch Processing",
                "description": "Load terrains from multiple projects",
                "code": """import arcpy
import os

# Directory containing HEC-RAS projects
project_dir = r"C:\\RAS_Projects"

# Find all project files
for root, dirs, files in os.walk(project_dir):
    for file in files:
        if file.endswith(".prj"):
            project_path = os.path.join(root, file)
            
            try:
                # Load primary terrain from each project
                arcpy.RASCommander.LoadRASTerrain(
                    in_ras_project=project_path,
                    import_all=False,
                    terrains_to_load=["Primary Terrain"]
                )
                print(f"Loaded terrain from: {file}")
            except Exception as e:
                print(f"Failed to load from {file}: {str(e)}")"""
            }
        ]
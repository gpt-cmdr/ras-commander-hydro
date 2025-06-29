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
    """
    def __init__(self):
        self.label = "Load HEC-RAS Terrain"
        self.description = "Loads terrain layers defined in a HEC-RAS project's .rasmap file into the current map."
        self.canRunInBackground = False
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
        params[1].value = False
        params[2].filter.type = "ValueList"
        
        return params

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
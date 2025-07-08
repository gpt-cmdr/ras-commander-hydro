# -*- coding: utf-8 -*-
#
# RAS Commander.pyt
#
# ArcGIS Python Toolbox for HEC-RAS HDF5 Data Integration
# ===================================================================================
#
# ESRI PARTNERSHIP ANNOUNCEMENT:
# ESRI has partnered with CLB Engineering Corporation's William "Bill" Katzenmeyer,
# Vice President and creator of the RAS Commander Open Source Python Library, to bring
# powerful HEC-RAS 6.x HDF5 data extraction capabilities directly into ArcGIS Pro.
#
# LAUNCHING AT ESRI USER CONFERENCE 2025
# This toolbox represents a groundbreaking application of LLM Forward engineering,
# developed in just over a month following an ASFPM brainstorming session in May 2025.
#
# DESCRIPTION:
# This toolbox provides comprehensive tools for loading and visualizing HEC-RAS 1D and 2D
# geometry, terrain, and results data from HDF5 files directly within ArcGIS Pro.
#
# KEY FEATURES:
# • Direct HDF5 data access without manual conversion
# • Support for 1D and 2D geometry elements
# • Pipe network extraction (storm/sewer systems)
# • Results visualization (Max WSE and Velocity)
# • Terrain import from RAS Mapper
# • Complete project organization
#
# ORIGIN AND ATTRIBUTION:
# This toolbox is a direct port of the HDF5 data extraction logic from the
# ras-commander library, adapted for the ArcGIS platform using CLB's innovative
# LLM Forward approach.
#
# RESOURCES AND LINKS:
# • RAS Commander Arc Hydro Tools: https://github.com/gpt-cmdr/ras-commander-hydro
# • RAS Commander Library: https://github.com/gpt-cmdr/ras-commander
# • CLB Engineering Corporation: https://clbengineering.com/
# • LLM Forward Approach: https://clbengineering.com/
# • Engineering with LLMs: https://engineeringwithllms.info/
#
# COMMUNITY DRIVEN:
# This is a community-driven effort. We're actively seeking feedback from:
# - Municipalities integrating HEC-RAS data into dashboards
# - Engineers communicating multi-hazard flood risk
# - GIS professionals preparing 2D model data
# - Researchers analyzing model results
#
# Share your ideas: https://github.com/gpt-cmdr/ras-commander-hydro/issues
#
# ===================================================================================

import sys
import os

# Add the Scripts directory to the Python path so we can import our modules
toolbox_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.join(os.path.dirname(toolbox_dir), 'Scripts', 'archydro')
if scripts_dir not in sys.path:
    sys.path.insert(0, scripts_dir)

# Import the tool classes from our modules
from rc_load_ras_terrain import LoadRASTerrain
from rc_load_hecras_2d_geometry import LoadHECRAS2DGeometry
from rc_load_hecras_2d_results import LoadHECRAS2DResults
from rc_load_hecras_1d_geometry import LoadHECRAS1DGeometry
from rc_organize_ras_project import OrganizeRASProject


class Toolbox(object):
    """
    ArcGIS Python Toolbox for loading HEC-RAS 1D and 2D geometry, terrain, and results layers.
    
    ESRI USER CONFERENCE 2025 LAUNCH EDITION
    
    Developed through ESRI's partnership with CLB Engineering Corporation's William "Bill" Katzenmeyer,
    this toolbox brings the power of the RAS Commander library directly into Arc Hydro Tools.
    """
    def __init__(self):
        # Core toolbox properties
        self.label = "RAS Commander Tools"
        self.alias = "RASCommander"
        self.description = """RAS Commander Arc Hydro Tools - Bringing HEC-RAS 6.x Direct Data Access to ArcGIS
        
        🎉 LAUNCHING AT ESRI USER CONFERENCE 2025 🎉
        
        This toolbox is the result of ESRI's partnership with CLB Engineering Corporation's 
        William "Bill" Katzenmeyer, creator of the RAS Commander Open Source Python Library.
        
        Using CLB's innovative LLM Forward approach, this toolbox was developed in just over 
        a month, demonstrating the transformative potential of AI-assisted development in 
        the water resources sector.
        
        KEY CAPABILITIES:
        • Direct HDF5 Import - No conversion needed
        • 1D and 2D Geometry - Complete model extraction
        • Pipe Networks - Storm/sewer system support
        • Results Analysis - Max WSE and velocity visualization
        • Terrain Integration - RAS Mapper VRT import
        • Project Organization - Batch processing tools
        
        COMMUNITY DRIVEN:
        We're seeking feedback to shape future development! Are you a:
        • Municipality looking to integrate HEC-RAS data into dashboards?
        • Engineer communicating multi-hazard flood risk?
        • GIS Professional preparing 2D model data?
        • Researcher analyzing model results?
        
        Share your use cases and feature requests!
        
        RESOURCES:
        • RAS Commander Arc Hydro Tools: https://github.com/gpt-cmdr/ras-commander-hydro
        • RAS Commander Library: https://github.com/gpt-cmdr/ras-commander
        • CLB Engineering: https://clbengineering.com/
        • Engineering with LLMs: https://engineeringwithllms.info/
        
        Developed by CLB Engineering in partnership with ESRI."""
        
        # Tool list
        self.tools = [LoadHECRAS1DGeometry, LoadHECRAS2DGeometry, LoadHECRAS2DResults, LoadRASTerrain, OrganizeRASProject]
        
        # Toolbox metadata
        self.author = "William 'Bill' Katzenmeyer, P.E. - CLB Engineering Corporation"
        self.credits = """ESRI PARTNERSHIP:
        Developed through ESRI's partnership with CLB Engineering Corporation
        
        DEVELOPMENT APPROACH:
        Created using CLB's LLM Forward methodology - Learn more at https://clbengineering.com/
        
        BASED ON:
        RAS Commander library: https://github.com/gpt-cmdr/ras-commander
        
        LEARN MORE:
        • About the author's work: https://engineeringwithllms.info/
        • CLB Engineering: https://clbengineering.com/"""
        
        self.version = "1.0.0 - Esri UC 2025 Launch Edition"
        self.homepage = "https://github.com/gpt-cmdr/ras-commander-hydro"
        
        # Keywords for searchability
        self.tags = ["HEC-RAS", "Arc Hydro", "HDF5", "River Analysis", "Hydraulic Modeling", 
                     "Flood Modeling", "2D Mesh", "RAS Mapper", "CLB Engineering", "LLM Forward",
                     "Esri UC 2025", "William Katzenmeyer", "Community Driven"]
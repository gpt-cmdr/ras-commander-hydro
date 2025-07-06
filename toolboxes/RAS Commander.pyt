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
# library's HDF handling classes (e.g., HdfMesh, HdfBndry).
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
    """
    def __init__(self):
        self.label = "Arc Hydro RAS-Commander Tools"
        self.alias = "RASCommander"
        self.description = "Tools for loading HEC-RAS 1D and 2D geometry, terrain, and results from HDF5 files. Sponsored by CLB Engineering (https://clbengineering.com/)."
        # List the tool classes
        self.tools = [LoadHECRAS1DGeometry, LoadHECRAS2DGeometry, LoadHECRAS2DResults, LoadRASTerrain, OrganizeRASProject]
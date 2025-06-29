# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is the Arc Hydro RAS Commander Tools repository - an ArcGIS Python Toolbox for HEC-RAS HDF5 data integration. The toolbox provides tools for loading and visualizing HEC-RAS 2D geometry and terrain data from HDF5 files directly within ArcGIS Pro.

**Sponsorship**: CLB Engineering (https://clbengineering.com/) in cooperation with ESRI

## Why Python Toolbox (.pyt)?

This project uses Python Toolbox format instead of binary toolbox format (.tbx) for several reasons:
- **System Toolbox Compatibility**: Better suited for system toolboxes in ArcGIS Pro
- **Version Control**: Text-based format ideal for Git (track changes, review diffs, collaborate)
- **No ArcGIS Pro Required**: Can be edited with any text editor
- **Dynamic Import**: Can import tool classes from separate Python modules
- **Cross-Platform**: Works identically across different OS and ArcGIS Pro installations

## Architecture

### Core Components

1. **Toolbox Structure** (refactored per Refactor Plan.txt):
   - Python Toolbox file: `toolboxes/RAS Commander.pyt`
   - Python scripts: `Scripts/ras_commander/` directory
     - `LoadRASTerrain.py` - Tool for loading HEC-RAS terrain layers
     - `LoadHECRAS6xHDFData.py` - Tool for loading HEC-RAS 6.x HDF data
     - `utils.py` - Shared utility functions
   - Original single-file toolbox: `ras-commander-hydro.pyt` (deleted - see git history)
   - Note: Production will also include an .atbx file that imports the .pyt classes

2. **Production Deployment Structure** (per Refactor Plan.txt):
   - Scripts: `C:\Program Files\ArcGIS\Pro\Resources\ArcToolBox\Scripts\ras_commander`
   - Toolbox: `C:\Program Files\ArcGIS\Pro\Resources\ArcToolBox\toolboxes\RAS Commander.pyt`
   - Layers: `C:\Program Files\ArcGIS\Pro\Resources\ArcToolBox\Templates\Layers\archydro\ras-commander`
   - Images: `C:\Program Files\ArcGIS\Pro\Resources\ArcToolBox\Images`
   - Schemas: `C:\Program Files\ArcGIS\Pro\Resources\ArcToolBox\Data\archydro\Ras2DTemplate.gdb`

3. **GitHub Repository Structure** (mirrors ArcToolbox folder):
   - `Scripts/ras_commander/` - Python tool implementations
   - `toolboxes/RAS Commander.pyt` - Main Python toolbox file
   - `Templates/Layers/archydro/ras-commander/` - Layer templates (.lyrx files)
   - `Images/` - Images (used sparingly, mainly for CLB Engineering logo)
   - `Data/archydro/Ras2DTemplate.gdb` - Geodatabase schema templates
   - `unittests/` - Unit tests for the Python scripts
   - `testdata/` - HDF test files for different scenarios

### Key Technical Details

- **HDF5 Processing**: Uses h5py and numpy for efficient HDF5 file reading
- **Projection Handling**: Extracts projection info from HDF attributes or associated .prj files
- **Geometry Creation**: Optimized polygon creation using numpy vectorization
- **ArcGIS Integration**: Uses arcpy for all GIS operations and messaging


## Commands

### Development Installation
```powershell
# Option 1: Run as administrator using batch file
install_toolbox_as_admin.bat

# Option 2: Run PowerShell script directly (requires admin)
.\install_toolbox.ps1

# The script will copy all components to the ArcGIS Pro installation directories
# When prompted, choose 'y' for development mode to create symlinks instead of copying
```

### Testing
```bash
# Run unit tests from the unittests folder
python -m pytest unittests/

# Test data files for different scenarios:
# - 1D Unsteady: testdata/BaldEagle.p01.hdf
# - 1D Steady: testdata/BaldEagle.p02.hdf
# - 2D Unsteady: testdata/BaldEagleDamBrk.p07.hdf
# - 2D with Pipes/Pumps: testdata/DavisStormSystem.p02.hdf
# - 2D with SWMM Import: testdata/BeaverLakeSWMMImpor.p01.hdf (contains pipe networks with field name issues)
```

## Important Notes

- **HEC-RAS Version**: The toolbox only supports HEC-RAS 6.x Model Series Results
- **Logging**: Helper functions use `arcpy.AddMessage/AddWarning` for logging since they don't have access to the tool's messages object
- **Utility Folders**: The `.ai_tools` folder contains LLM knowledge base generation utilities - ignore it for most purposes
- **Layer Templates**: The repository includes pre-configured layer templates (.lyrx files) for visualizing different HEC-RAS data types:
  - 2D Breaklines
  - Maximum Velocity at Cell Faces
  - Maximum Water Surface Elevation at Cell Centers
  - Mesh Area Perimeters
  - Mesh Cell Centers, Faces, and Polygons
  - Pipe Conduits and Nodes (for storm/sewer network visualization)

## HDF Data Field Handling

### Known Issues and Solutions

1. **Field Name Conflicts with ArcGIS System Fields**
   - Fields named "Shape", "OBJECTID", "SHAPE_LENGTH", etc. conflict with ArcGIS system fields
   - Solution: The tool automatically renames these fields (e.g., "Shape" → "Shape_Type")

2. **Special Characters in Field Names**
   - HEC-RAS field names may contain characters that ArcGIS doesn't support (apostrophes, colons, parentheses)
   - Solution: The tool automatically cleans field names by replacing special characters with underscores and removing apostrophes

3. **Field Name Typos in HDF Files**
   - Some HDF files contain typos (e.g., "Condtui Connections" instead of "Conduit Connections")
   - Solution: The tool automatically corrects known typos during field name processing

4. **NaN Values in Numeric Fields**
   - HDF files may contain NaN (Not a Number) values that cause type compatibility errors
   - Solution: The tool converts NaN values to NULL for proper handling in ArcGIS

5. **Boolean Fields**
   - HDF boolean fields need special handling for ArcGIS compatibility
   - Solution: Boolean values are converted to SHORT fields (0/1 values)

### Field Type Detection
The tool automatically detects field types based on data values:
- **SHORT**: Boolean values
- **LONG**: Integer values
- **DOUBLE**: Floating-point values (including fields with all NaN values that sound numeric)
- **TEXT**: String values (with calculated field lengths based on actual data)

### Debugging Field Issues
When encountering field-related errors:
1. Check the tool messages for DEBUG output showing original vs. cleaned field names
2. Look for field mapping errors that show expected vs. available fields
3. The tool provides detailed error messages when field type conversions fail

## Troubleshooting

### Common Issues
1. **"Script cannot be loaded"**: PowerShell execution policy issue
   - Use the batch file or run: `Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser`

2. **"Access denied"**: Ensure running as administrator and ArcGIS Pro is closed

3. **"ArcGIS Pro not found"**: Verify ArcGIS Pro installation in standard location

4. **"Shape already exists" error**: Field name conflicts with ArcGIS system fields
   - The tool should automatically handle this, but check DEBUG messages for field renaming

5. **"Cannot find field" error**: Field name mismatch between HDF data and ArcGIS
   - Check tool messages for field name transformations
   - Look for typos in HDF field names that need correction

6. **"Value type is incompatible with field type" error**: NaN or type conversion issues
   - Usually caused by NaN values in numeric fields
   - The tool automatically converts NaN to NULL

### Verifying Installation
1. Open ArcGIS Pro
2. Go to Catalog pane → Toolboxes → System Toolboxes
3. Look for "RAS Commander"
4. Tools available: Load HEC-RAS Terrain, Load HEC-RAS 6.x HDF Data
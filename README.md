# RAS Commander Arc Hydro Tools

<div align="center">
  <img src="Images/ras-commander-archydro.svg" alt="RAS Commander Arc Hydro Tools" width="600">
  
  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
  [![ArcGIS Pro](https://img.shields.io/badge/ArcGIS%20Pro-2.8%2B-blue)](https://www.esri.com/en-us/arcgis/products/arcgis-pro/overview)
  [![HEC-RAS](https://img.shields.io/badge/HEC--RAS-6.x-green)](https://www.hec.usace.army.mil/software/hec-ras/)
  
  **Bringing HEC-RAS 6.x HDF5 Data to ArcGIS Pro**
  
  [CLB Engineering](https://clbengineering.com/) | [RAS Commander Library](https://github.com/gpt-cmdr/ras-commander) | [Arc Hydro](https://www.esri.com/en-us/industries/water-resources/arc-hydro)
</div>

---

## Overview

**RAS Commander Arc Hydro Tools** is an ArcGIS Python Toolbox that brings powerful HEC-RAS 6.x HDF5 data extraction capabilities directly into ArcGIS Pro. This toolbox enables hydraulic engineers and GIS professionals to seamlessly load and visualize HEC-RAS 1D and 2D geometry, terrain, and results data without manual conversion steps.

### Key Features

- 📊 **Direct HDF5 Import** - Load HEC-RAS data directly from geometry (g*.hdf) and plan (p*.hdf) files
- 🗺️ **1D Geometry Support** - Extract cross sections, river centerlines, bank lines, and structures
- 🌊 **2D Geometry Support** - Import mesh elements, breaklines, boundary conditions, and cell polygons
- 🏗️ **Pipe Networks** - Full support for storm/sewer pipe networks including SWMM imports
- 📈 **Results Visualization** - Display maximum WSE and velocity results with time of occurrence
- ⛰️ **Terrain Loading** - Import HEC-RAS terrain layers from RASMapper VRT files
- 🗂️ **Project Organization** - Batch process entire HEC-RAS projects into organized geodatabases

<div align="center">
  <a href="https://clbengineering.com/">
    <img src="Images/CLBEngineeringMainLogo.png" alt="CLB Engineering" width="300">
  </a>
  
  *Sponsored by [CLB Engineering](https://clbengineering.com/) in cooperation with ESRI*
</div>

---

## Installation

### Primary Method: Arc Hydro Tools Installation

The RAS Commander toolbox is included as part of the Arc Hydro Tools distribution. This is the recommended installation method for most users.

1. Install Arc Hydro Tools following the standard installation process
2. The RAS Commander toolbox will be available under:
   ```
   Toolboxes → Arc Hydro Tools → RAS Commander
   ```

### Development Installation

For developers and users who want to extend or customize the tools:

1. **Fork and Clone the Repository**
   ```bash
   git clone https://github.com/your-fork/ras-commander-hydro.git
   cd ras-commander-hydro
   ```

2. **Option A: Add Toolbox in ArcGIS Pro**
   - Open ArcGIS Pro
   - In the Catalog pane, right-click on Toolboxes
   - Select "Add Toolbox"
   - Navigate to `toolboxes/RAS Commander.pyt`

3. **Option B: Install Permanently (Development Only)**
   ```bash
   # Run as Administrator
   install_toolbox_as_admin.bat
   ```
   
   ⚠️ **Warning**: This installation will be overwritten if Arc Hydro Tools is installed or updated. This method is for development only.

---

## Tools Overview

### 🔧 Load HEC-RAS 1D Geometry Layers
Extract 1D hydraulic model elements including:
- Cross sections with station-elevation data
- River and reach centerlines
- Left and right bank lines
- Edge lines for terrain processing
- Bridges, culverts, weirs, and other structures

### 🌐 Load HEC-RAS 2D Geometry Layers
Import 2D model components including:
- 2D flow area perimeters
- Mesh cell centers, faces, and polygons
- Breaklines with cell spacing attributes
- External and internal boundary conditions
- Pipe conduits and junction nodes (storm/sewer networks)

### 📊 Load HEC-RAS 2D Results Summary Layers
Visualize simulation results including:
- Maximum water surface elevation at cell centers
- Maximum velocity at cell faces
- Time of maximum occurrence for all results

### ⛰️ Load HEC-RAS Terrain
Import terrain layers from RASMapper:
- Loads VRT (Virtual Raster) files with HEC-RAS priority
- Supports multiple terrain layers per project

⚠️ **Important**: Only base terrain rasters are loaded. Vector terrain modifications (breaklines, high ground, etc.) made in RAS Mapper are not included.

### 🗂️ Organize HEC-RAS Project
Comprehensive project processing tool that:
- Processes all plan files in a project directory
- Creates organized geodatabase structure
- Extracts all available geometry and results
- Groups outputs by project and plan number
- Automatically adds results to current map

---

## Development Narrative: AI-Powered Tool Creation

### The Power of LLM-Assisted Development

The RAS Commander Arc Hydro Tools represent a unique achievement in AI-assisted software development. This toolbox was created through an innovative process that demonstrates the transformative potential of Large Language Models (LLMs) in software engineering.

#### Foundation: The RAS Commander Library

The journey began with the [ras-commander library](https://github.com/gpt-cmdr/ras-commander), an LLM-driven automation library built as a modern Python complement to the classic HECRASController. This library included several HDF5-to-GeoPandas conversion functions that served as the foundation for our ArcGIS tools.

#### Transfer Learning Between Frameworks

What makes this development process remarkable is the LLM's ability to perform **transfer learning** between different frameworks and constraints:

1. **Framework Translation**: The LLM successfully translated GeoPandas-based operations to ArcPy equivalents while maintaining functional integrity
2. **Constraint Adaptation**: The code was adapted to work within ArcGIS Pro's Python environment without adding any dependencies beyond the standard ArcPy installation
3. **Metadata Preservation**: Where possible, the tools preserve HEC-RAS HDF5 data structures and naming conventions for consistency

#### Self-Contained Architecture

A key achievement was creating completely self-contained tools that:
- Rely only on standard ArcPy and h5py (included with ArcGIS Pro)
- Require no additional package installations
- Maintain high performance through optimized numpy operations
- Handle complex geometries and large datasets efficiently

### Technical Innovations

#### Dynamic Field Mapping

One of the most interesting challenges involved handling pipe network data from SWMM imports. HEC-RAS stores these with column names that sometimes include typos or unconventional formatting. The tools dynamically:
- Detect and correct field name typos (e.g., "Condtui_Connections" → "Conduit_Connections")
- Handle reserved ArcGIS field names (e.g., "Shape" → "Shape_Type")
- Preserve all original HDF5 attributes while ensuring ArcGIS compatibility

#### Intelligent Data Extraction

The tools implement smart extraction strategies:
- Pre-compute mesh faces for efficient polygon construction
- Use vectorized numpy operations for coordinate transformations
- Implement optimized batch insertion for large datasets
- Cache HDF5 metadata to minimize file access

### The LLM Advantage

This project demonstrates several unique advantages of LLM-assisted development:

1. **Rapid Prototyping**: Complex tool logic was developed and refined through iterative conversations
2. **Cross-Domain Knowledge**: The LLM could apply hydraulic modeling concepts while respecting GIS constraints
3. **Best Practices Integration**: Modern Python patterns and ArcPy best practices were automatically incorporated
4. **Documentation Generation**: Comprehensive help documentation and metadata were created alongside the code

The result is a professional-grade toolbox that bridges two specialized domains—hydraulic modeling and GIS—through the power of AI-assisted development.

---

## Usage Examples

### Basic 1D Geometry Import
```python
import arcpy

# Load cross sections and river centerlines
arcpy.RASCommander.LoadHECRAS1DGeometry(
    input_hdf=r"C:\Models\MyProject.g01.hdf",
    geometry_elements=["Cross Sections", "River Centerlines"],
    output_cross_sections=r"C:\Output\Output.gdb\CrossSections",
    output_centerlines=r"C:\Output\Output.gdb\RiverCenterlines"
)
```

### Organize Entire Project
```python
# Process all plan files in a project
arcpy.RASCommander.OrganizeRASProject(
    input_path=r"C:\Models\MyProject",
    output_gdb=r"C:\Output\MyProject_Organized.gdb",
    include_cell_polygons=True,
    extract_all_results=True
)
```

---

## Requirements

- **ArcGIS Pro**: Version 2.8 or higher
- **HEC-RAS**: Version 6.x model files
- **Python**: Uses ArcGIS Pro's built-in Python environment
- **Dependencies**: None beyond standard ArcPy (h5py and numpy are included)

---

## Project Structure

```
ras-commander-hydro/
├── Images/                          # Branding and documentation images
├── Scripts/
│   └── archydro/                   # Tool implementation modules
│       ├── rc_load_hecras_1d_geometry.py
│       ├── rc_load_hecras_2d_geometry.py
│       ├── rc_load_hecras_2d_results.py
│       ├── rc_load_ras_terrain.py
│       ├── rc_organize_ras_project.py
│       └── rc_utils.py             # Shared utilities
├── toolboxes/
│   ├── RAS Commander.pyt           # Main Python toolbox
│   └── *.pyt.xml                   # Tool metadata files
├── Doc/
│   └── RASCommander_Help.html      # Integrated help documentation
├── install_toolbox.ps1             # PowerShell installation script
├── install_toolbox_as_admin.bat    # Admin installation batch file
└── README.md                       # This file
```

---

## Contributing

We welcome contributions! The RAS Commander Arc Hydro Tools are part of a community-driven effort to improve hydraulic modeling workflows.

### How to Contribute

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Development Guidelines

- Follow the existing code style and naming conventions
- Maintain compatibility with Arc Hydro Tools standards
- Use the `rc_` prefix for all Python modules
- Test with various HEC-RAS model types (1D, 2D, combined)
- Update documentation for new features

---

## Future Roadmap

### Planned Features

- **1D Results**: Water surface profiles, velocities, and flow data
- **Time Series Results**: Full temporal data extraction and animation
- **Floodplain Mapping**: Automated flood extent and depth grid generation
- **Model Synchronization**: Two-way data exchange with RAS Mapper
- **Cloud Integration**: Support for cloud-hosted HEC-RAS models
- **Performance Monitoring**: Model runtime statistics and diagnostics

### Community Requests

We're actively seeking feedback from:
- 🏛️ **Municipalities** looking to integrate HEC-RAS data into dashboards
- 👷 **Engineers** communicating multi-hazard flood risk
- 🗺️ **GIS Professionals** preparing 2D model data
- 🔬 **Researchers** analyzing model results

[Share your ideas and use cases!](https://github.com/gpt-cmdr/ras-commander-hydro/issues)

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **[CLB Engineering](https://clbengineering.com/)** - Project sponsor and hydraulic engineering expertise
- **[ESRI](https://www.esri.com/)** - Arc Hydro Tools integration and GIS platform
- **[USACE HEC](https://www.hec.usace.army.mil/)** - HEC-RAS software and HDF5 format documentation
- **[ras-commander](https://github.com/gpt-cmdr/ras-commander)** - Original Python library and HDF5 extraction logic

## Trademarks

- **HEC-RAS™** is a trademark of the U.S. Army Corps of Engineers (USACE) Hydrologic Engineering Center (HEC)
- **ARC HYDRO** is a trademark of Environmental Systems Research Institute (ESRI)
- **RAS Commander™** is a trademark of CLB Engineering Corporation

*RAS Commander Arc Hydro Tools is an independent project and is not affiliated with, endorsed by, or sponsored by USACE or HEC.*

---

<div align="center">
  <img src="Images/ras-commander_logo.svg" alt="RAS Commander" width="150">
  
  **Transform Your HEC-RAS Workflow Today**
  
  [Get Started](https://github.com/gpt-cmdr/ras-commander-hydro) | [Documentation](Doc/RASCommander_Help.html) | [Report Issues](https://github.com/gpt-cmdr/ras-commander-hydro/issues)
</div>
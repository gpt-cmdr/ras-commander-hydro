# RAS Commander for ArcGIS

![CLB Engineering Logo](https://clbengineering.com/wp-content/uploads/2023/11/CLB-Horizontal-Logo-Color-e1700244944888.png)

## Overview

**RAS Commander for ArcGIS** is a powerful Python Toolbox (`.pyt`) designed to seamlessly integrate HEC-RAS data into your ArcGIS Pro projects. It provides a suite of tools to extract and load 2D geometry layers directly from HEC-RAS model HDF5 files (`.g*.hdf`, `.p*.hdf`) into geodatabase feature classes.

This toolbox eliminates the need for manual data conversion or intermediate file formats, creating a direct and efficient workflow for hydraulic modelers and GIS professionals.

### Sponsorship

The development of this ArcGIS toolbox was generously sponsored by **[CLB Engineering](https://clbengineering.com/)** in cooperation with **ESRI**. Their support makes this open-source contribution to the community possible.

### Origin and Attribution

This toolbox is a direct port of the HDF5 data extraction logic from the **[ras-commander library](https://github.com/gpt-cmdr/ras-commander)**. All core HDF5 reading logic is derived from the library's HDF handling classes (e.g., `HdfMesh`, `HdfBndry`). We extend our thanks to the original developers for their foundational work.

## Features

- **Direct HDF5 Reading:** No need to export HEC-RAS data to other formats.
- **Self-Contained:** The `.pyt` file includes all necessary logic and has no external dependencies beyond a standard ArcGIS Pro installation.
- **Automatic CRS Detection:** The tools attempt to automatically detect the Coordinate Reference System (CRS) from the HEC-RAS project files.
- **CRS Override:** Provides an option to manually specify a CRS if automatic detection fails.
- **Robust and Efficient:** Uses `arcpy.da.InsertCursor` for efficient writing of features, capable of handling very large datasets.

## Included Tools

The toolbox provides the following tools:

- **Load 2D Breaklines:** Extracts mesh breaklines as a polyline feature class.
- **Load 2D Mesh Area Perimeters:** Extracts the outer boundary of each 2D flow area as a polygon feature class.
- **Load 2D Mesh Cell Centers:** Extracts the center point of every cell in a 2D mesh as a point feature class.
- **Load 2D Mesh Cell Faces:** Extracts the faces (edges) between mesh cells as a polyline feature class.
- **Load 2D Mesh Cells as Polygons:** Reconstructs and loads the full polygon geometry for every cell in a 2D mesh.
- **About RAS Commander:** Displays information about the toolbox, its origin, and sponsorship.

## Installation and Usage

1.  **Download:** Obtain the `ras-commander-hydro.pyt` file.
2.  **Add to ArcGIS Pro:**
    -   Open ArcGIS Pro and go to the **Catalog** pane.
    -   Right-click on **Toolboxes** and select **Add Toolbox**.
    -   Navigate to and select the `ras-commander-hydro.pyt` file.
3.  **Run a Tool:**
    -   Expand the "RAS Commander (by CLB Engineering)" toolbox in the Catalog pane.
    -   Double-click on any tool to open its dialog.
    -   Fill in the required parameters (e.g., input HDF file, output feature class location).
    -   Click **Run**.

## Requirements

-   **ArcGIS Pro:** Version 2.5 or newer.
-   **HEC-RAS HDF Files:** Geometry (`.g*.hdf`) or Plan (`.p*.hdf`) files from HEC-RAS 6.0 or newer.

## License

This software is licensed under the MIT License. See the `LICENSE` file for more details.

## Contributing

Contributions are welcome! If you have suggestions for improvements or new features, please open an issue or submit a pull request on the [GitHub repository](https://github.com/gpt-cmdr/ras-commander-hydro).
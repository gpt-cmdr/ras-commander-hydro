# RAS Commander Hydro

!!! note "A focused tool, not the full library"
    RAS Commander Hydro provides a **focused subset** of helper tools. For the most powerful and complete way to automate HEC-RAS, use the [**ras-commander** library](https://rascommander.info/ras/) directly.

**RAS Commander Arc Hydro Tools** brings HEC-RAS 6.x direct data access into ArcGIS Pro. It
lets hydraulic engineers and GIS professionals extract and visualize HEC-RAS 1D and 2D
geometry, terrain, and results data without manual conversion steps. The toolbox is free and
open source, and is built on the [ras-commander library](https://github.com/gpt-cmdr/ras-commander).

The toolbox was developed by [CLB Engineering](https://clbengineering.com/) in partnership with
ESRI, integrating into the [Arc Hydro Tools](https://www.esri.com/en-us/industries/water-resources/arc-hydro)
framework using CLB's LLM Forward approach.

## Key Features

- **Direct HDF5 Import** — Load HEC-RAS data directly from geometry (`g*.hdf`) and plan
  (`p*.hdf`) files.
- **1D Geometry Support** — Extract cross sections, river centerlines, bank lines, and
  structures.
- **2D Geometry Support** — Import mesh elements, breaklines, boundary conditions, and cell
  polygons.
- **Pipe Networks** — Full support for storm/sewer pipe networks including SWMM imports.
- **Results Visualization** — Display maximum WSE and velocity results with time of occurrence.
- **Terrain Loading** — Import HEC-RAS terrain layers from RAS Mapper VRT files.
- **Project Organization** — Batch process entire HEC-RAS projects into organized geodatabases.

## Who It's For

- **Municipalities** integrating HEC-RAS data into dashboards.
- **Engineers** communicating multi-hazard flood risk.
- **GIS professionals** preparing 2D model data.
- **Researchers** analyzing model results.

## Where to Go Next

- [Installation](installation.md) — install via Arc Hydro Tools or set up a development copy.
- [Features](features.md) — a reference of the available tools.

## Resources

- [RAS Commander Arc Hydro Tools (GitHub)](https://github.com/gpt-cmdr/ras-commander-hydro)
- [RAS Commander Library](https://github.com/gpt-cmdr/ras-commander)
- [CLB Engineering Corporation](https://clbengineering.com/)
- [Engineering with LLMs](https://engineeringwithllms.info/)

Licensed under the MIT License.

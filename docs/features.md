# Features

RAS Commander Arc Hydro Tools provides a focused set of tools for importing HEC-RAS 6.x data
into ArcGIS Pro. Each tool reads directly from HEC-RAS geometry (`g*.hdf`) and plan (`p*.hdf`)
files.

## Load HEC-RAS 1D Geometry Layers

Extract comprehensive 1D hydraulic model elements — cross sections, river centerlines, bank
lines, and structures — for report figures and analysis.

<p align="center">
  <img src="../Images/docs/Load1DGeometry.png" alt="Load 1D Geometry" width="35%">
</p>

## Load HEC-RAS 2D Geometry Layers

Import complete 2D model components including mesh elements, breaklines, boundary conditions,
and mesh cells as polygons for advanced spatial analysis. Includes full support for storm/sewer
pipe networks (including SWMM imports).

<p align="center">
  <img src="../Images/docs/Load2DGeometry.png" alt="Load 2D Geometry" width="35%">
</p>

## Load HEC-RAS 2D Results Summary Layers

Visualize maximum water-surface elevation and velocity results, with time of occurrence.

<p align="center">
  <img src="../Images/docs/Load2DSummaryResults.png" alt="Load 2D Results Summary" width="35%">
</p>

## Load HEC-RAS Terrain

Import terrain layers from RAS Mapper VRT files with proper georeferencing.

<p align="center">
  <img src="../Images/docs/LoadRASTerrain.png" alt="Load RAS Terrain" width="35%">
</p>

## Organize HEC-RAS Project

A batch-processing tool that organizes an entire HEC-RAS project into an organized geodatabase.

<p align="center">
  <img src="../Images/docs/OrganizeRASProject.png" alt="Organize HEC-RAS Project" width="35%">
</p>

## Example: New Orleans 2D Model

The images below show a 2D HEC-RAS model of the New Orleans metro storm-water system, complete
with imported pipe networks, mesh polygons, and a maximum WSEL raster generated directly inside
ArcGIS Pro.

<p align="center">
  <img src="../Images/docs/rc_neworleanspipes.png" alt="New Orleans Imported Pipe Networks" width="50%">
</p>

<p align="center">
  <img src="../Images/docs/rc_neworleanspipes_results.png" alt="New Orleans Maximum WSEL" width="50%">
</p>

## Roadmap

Initial release capabilities:

- 1D and 2D geometry extraction (including pipe networks).
- Max WSE and velocity as 2D mesh results.
- Terrain import for inundation mapping.
- Support for HEC-RAS 2D models.
- Organize entire projects as geodatabases.

Coming soon:

- Improved schemas and layer styling.
- 1D results and full time series.
- Fluvial/pluvial delineation.
- Land use layer integration.
- Sync changes back to HEC-RAS.

[View the full roadmap and vote on features](https://github.com/gpt-cmdr/ras-commander-hydro/issues).

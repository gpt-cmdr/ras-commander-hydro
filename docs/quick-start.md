# Quick Start

This workflow produces a reviewable geodatabase while protecting the original HEC-RAS project.

## 1. Open a map and identify the model files

Work from a copy or read-only source model when practical. In the HEC-RAS project folder, identify:

| File | Use |
| --- | --- |
| `Project.g01.hdf` | Geometry HDF; suitable for 1D or 2D geometry tools. |
| `Project.p01.hdf` | Plan HDF; suitable for geometry and, after a successful compute, summary results. |
| `Project.prj` | Project file used by the terrain loader. |
| `Project.rasmap` | RAS Mapper configuration read automatically by the terrain loader. |
| `Terrain\...\Terrain.vrt` | Base terrain virtual raster referenced from the `.rasmap` file. |

Plan and geometry numbers vary. A computed plan HDF is normally the clearest single input because
it ties geometry and results to one plan.

## 2. Choose an output strategy

For a durable deliverable, create a file geodatabase on a local drive. The individual geometry and
results tools default to creating `Project.pXX.gdb` beside the selected HDF. Within it, outputs are
placed in a feature dataset named `Project_PlanXX`.

Use in-memory outputs only for quick inspection. They are temporary and are a poor choice for large
meshes or a deliverable.

## 3. Confirm the coordinate system

The geometry and results tools first try to read the HEC-RAS projection from the HDF/project
context. **Override CRS** is used only when automatic detection fails; it does not reproject data
that already has a detected coordinate system.

!!! warning
    Do not choose a CRS merely to make a layer draw near the project. Confirm it from the HEC-RAS
    project, RAS Mapper, survey/control data, or other authoritative project documentation.

## 4. Run the appropriate tool

- For one layer type, use an [individual loader](features.md).
- For the available geometry and summary results across a project folder, use
  [Organize HEC-RAS Project](tools/organize-project.md).

Start with mesh-area perimeters or center points before requesting cell polygons from a large 2D
model. Cell-polygon construction can be much slower and more storage-intensive.

## 5. Review the messages and outputs

After the tool finishes:

1. Open **Geoprocessing History** and read all warnings.
2. Confirm that the reported feature counts are plausible.
3. Check the output coordinate system and overlay against trusted project data.
4. Confirm that the plan number and HDF path match the intended simulation.
5. For results, confirm units and compare several values with RAS Mapper.
6. Save the ArcGIS Pro project if you want to retain map organization and symbology.

Use the full [validation checklist](reference/validation.md) before sharing outputs.

## Recommended first runs

### Focused 2D review

1. Run **Load HEC-RAS 2D Geometry Layers** for **Mesh Area Perimeters** and **Mesh Cell Centers**.
2. Run **Load HEC-RAS 2D Results Summary Layers** for **Max WSE at Cell Centers**.
3. Compare the output extent and sample values with RAS Mapper.
4. Add faces, polygons, breaklines, or pipe layers only as needed.

### Complete project inventory

1. Run **Organize HEC-RAS Project** on the project folder.
2. Keep 1D geometry, 2D geometry, and 2D results enabled.
3. Leave mesh-cell polygons disabled for the first run.
4. Review each plan group and rerun selectively if a plan needs deeper inspection.

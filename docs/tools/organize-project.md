# Organize HEC-RAS Project

Use this batch tool to extract available 1D geometry, 2D geometry, pipe networks, and 2D summary
results into one file geodatabase.

## Input behavior

The input can be a single HDF or a project folder.

- For a folder, the tool searches the folder's top level for `*.pNN.hdf` files.
- If no plan HDF files are found, it searches the same level for `*.gNN.hdf` files.
- The search is not recursive.
- For a single HDF, only that file is processed.

Place or select the intended plan HDFs directly rather than relying on copied results in nested
folders.

## Parameters

| Parameter | Default | Behavior |
| --- | --- | --- |
| **HEC-RAS Project Directory or HDF File** | Required | Project folder or one plan/geometry HDF. |
| **Output Geodatabase** | Derived from input | Required file geodatabase. A folder input defaults to `FolderName_Organized.gdb`; a file input defaults to `FileName_Organized.gdb`. |
| **Override CRS** | Empty | Used for files whose CRS cannot be determined automatically. |
| **Include 1D Geometry** | On | Extracts detected cross sections, centerlines, bank/edge lines, and structures. |
| **Include 2D Geometry** | On | Extracts detected mesh, boundary, breakline, and pipe elements. |
| **Include 2D Results Summary Layers** | On | Extracts maximum WSE and maximum face velocity when supported results exist. |
| **Include Mesh Cell Polygons** | Off | Adds reconstructed cell polygons; can substantially increase time and storage. |
| **Load Results to Map** | On | Adds created feature classes to the active map, grouped by plan. Disable for noninteractive or output-only runs. |

## Procedure

1. Open **Organize HEC-RAS Project**.
2. Select one plan HDF for a controlled run, or select the project folder for all top-level plans.
3. Review the automatically proposed output geodatabase.
4. Select the geometry/results categories required for the delivery.
5. Leave mesh-cell polygons disabled unless they are specifically needed.
6. Disable **Load Results to Map** when running without an active map or when the geodatabase alone
   is the deliverable.
7. Run the tool and review the availability report for every HDF.
8. Review each plan feature dataset and complete the [validation checklist](../reference/validation.md).

## Output organization

Each plan is written to a feature dataset based on the project and plan number, for example:

```text
Model_Organized.gdb
└── Model_Plan01
    ├── Model_Plan_01_CrossSections
    ├── Model_Plan_01_RiverCenterlines
    ├── Model_Plan_01_MeshPerimeters
    ├── Model_Plan_01_MeshCellCenters
    ├── Model_Plan_01_MaxWSE
    └── Model_Plan_01_MaxVelocity
```

Only available and selected data are written. A missing results group does not prevent available
geometry from being processed; the tool reports a warning and continues.

## Performance guidance

- Run one plan first to establish expected time, layer counts, and storage.
- Keep outputs on a local drive during extraction.
- Leave cell polygons off for large meshes; add them in a targeted second run if needed.
- Disable **Load Results to Map** for large batch runs to avoid the cost of adding every layer.
- Close unnecessary ArcGIS Pro views and confirm adequate disk space.

## Limitations

- Folder discovery is top-level and uses two-digit `pNN`/`gNN` naming.
- When any plan HDFs are present, geometry HDFs are not separately processed.
- The batch tool extracts only the five individual tools' supported layers.
- It does not copy source HEC-RAS files or create a self-contained HEC-RAS project.
- It does not execute plans, extract full time series, or write changes back to HEC-RAS.

See [Outputs and naming](../reference/outputs.md) for fields and naming details.

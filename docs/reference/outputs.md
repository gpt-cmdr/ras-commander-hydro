# Outputs and Naming

This reference describes the GIS products written by the current toolbox. The source HEC-RAS files
are opened for reading and are not modified.

## Individual-tool geodatabases

The 1D geometry, 2D geometry, and 2D results tools default to **Create New Geodatabase = True**.
For a plan HDF named `Model.p01.hdf`, the default organization is:

```text
Model.p01.gdb
└── Model_Plan01
    ├── CrossSections
    ├── MeshPerimeters
    ├── MaxWSE
    └── ...
```

The existing named geodatabase is reused if it already exists. Use a new output name when you need
to preserve an earlier extraction unchanged.

If geodatabase output is disabled, individual output parameters default to ArcGIS Pro's in-memory
workspace. In-memory layers are temporary and should not be treated as a deliverable.

## Batch-tool naming

**Organize HEC-RAS Project** prefixes feature-class names with project and plan identifiers, for
example `Model_Plan_01_MeshCellCenters`. This makes layers from multiple plans distinguishable in
one geodatabase and supports plan group layers in the map.

Plan parsing expects a base name ending in `.pNN`. Other names fall back to plan `00`; use standard
HEC-RAS names when predictable plan identification matters.

## 1D fields

| Feature class | Geometry | Fields written by the tool |
| --- | --- | --- |
| `CrossSections` | Polyline | `xs_id`, `River`, `Reach`, `RS`, `LeftBank`, `RightBank`, `LenLeft`, `LenChannel`, `LenRight`, `StationElevation`, `ManningsN` |
| `RiverCenterlines` | Polyline | `river_id`, `RiverName`, `ReachName`, `USType`, `DSType` |
| `BankLines` | Polyline | `bank_id`, `BankSide`, `Length` |
| `EdgeLines` | Polyline | `edge_id`, `EdgeType`, `Length` |
| `Structures1D` | Polyline | `struct_id`, `Type`, `River`, `Reach`, `RS`, `Description` |

ArcGIS adds standard fields such as `OBJECTID`, `Shape`, and geometry measurements as applicable.

## 2D fields

| Feature class | Geometry | Fields written by the tool |
| --- | --- | --- |
| `Breaklines2D` | Polyline | `bl_id`, `Name`, `CellSpaceNear`, `CellSpaceFar`, `NearRepeats`, `ProtectRadius` |
| `BCLines2D` | Polyline | `bc_id`, `Name`, `Type` |
| `MeshPerimeters` | Polygon | `mesh_name` |
| `MeshCellCenters` | Point | `mesh_name`, `cell_id` |
| `MeshCellFaces` | Polyline | `mesh_name`, `face_id` |
| `MeshCellPolygons` | Polygon | `mesh_name`, `cell_id` |
| `PipeConduits` | Polyline | Dynamic fields from available HDF conduit attributes |
| `PipeNodes` | Point | Dynamic fields from available HDF node attributes |
| `PipeNetworks` | Polygon | Dynamic fields from available HDF network attributes |

Dynamic pipe field names are normalized for ArcGIS. Spaces and punctuation may become underscores,
and names that conflict with ArcGIS fields may be renamed.

## Results fields

| Feature class | Geometry | Fields written by the tool |
| --- | --- | --- |
| `MaxWSE` | Point | `mesh_name`, `cell_id`, `max_wse`, `max_wse_time` |
| `MaxVelocity` | Point | `mesh_name`, `face_id`, `max_vel`, `time_of_max` |

Result dates are calculated from the HDF simulation start time plus the stored time-of-maximum
offset. Result values remain in source-model units.

## Coordinate systems

The toolbox assigns the CRS detected from the HEC-RAS project context. If none can be found, the
user must provide **Override CRS**. The override supplies a spatial reference; it is not a coordinate
transformation. Confirm alignment before analysis.

## Missing data

HEC-RAS HDF content varies by model type, version, and compute state. The toolbox writes an output
only when the required HDF group and datasets exist. A selected but unavailable element can produce
a message and no feature class or an empty feature class; review geoprocessing messages and counts.

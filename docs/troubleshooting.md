# Troubleshooting

Start with **Geoprocessing History** in ArcGIS Pro. The toolbox reports detected CRSs, HDF content,
skipped elements, missing files, and feature counts there.

## The toolbox does not appear

- Confirm that you installed the Arc Hydro build compatible with your ArcGIS Pro version.
- Restart ArcGIS Pro after installing Arc Hydro.
- For a source install, add `toolboxes/RAS-Commander.pyt`, not an individual script.
- Keep `toolboxes`, `Scripts/archydro`, and their relative repository layout intact.

## The tool reports that the CRS cannot be determined

1. Confirm the project projection in RAS Mapper and project documentation.
2. Check that associated project/projection files remain beside the HDF.
3. Set **Override CRS** to the confirmed spatial reference.
4. After extraction, overlay the result against trusted control data.

**Override CRS** assigns the coordinate system. It does not transform incorrectly projected
coordinates.

## The tool completes but an output is empty

- Check whether the selected element exists in the chosen geometry/plan.
- Try the plan HDF associated with the intended geometry.
- Review messages for missing HDF groups or incomplete datasets.
- Confirm that the output parameter was enabled and has a valid path.
- For pipe networks, confirm that the model actually contains the supported pipe HDF groups.

## No 2D results are found

- Use a computed `p*.hdf`, not a `g*.hdf`.
- Open the plan in HEC-RAS/RAS Mapper and confirm that 2D summary output exists.
- Confirm that the plan compute finished successfully and that the HDF is not stale or incomplete.
- Verify that the intended 2D flow area contains maximum water-surface or face-velocity summary datasets.

## The terrain list is empty

- Confirm that `Project.rasmap` exists beside `Project.prj` with the same base name.
- Open and save the project in RAS Mapper so terrain definitions are present.
- Confirm that the `.rasmap` contains terrain layer entries.
- Repair broken terrain paths in HEC-RAS/RAS Mapper rather than editing the GIS output.

## A terrain VRT is missing

The tool derives a `.vrt` path from the terrain path stored in `.rasmap`. Confirm that the matching
VRT exists at the reported location. If it does not, rebuild or repair the terrain through the
project's RAS Mapper workflow.

## The terrain differs from RAS Mapper

The tool loads only the base VRT. RAS Mapper vector terrain modifications are not applied. Confirm
whether the apparent difference is caused by high-ground, channel, levee, or other terrain
modifications in the HEC-RAS project.

## Cell polygons are slow or use too much storage

- Cancel only if safe for the current ArcGIS operation, then rerun without cell polygons.
- Use mesh perimeters, centers, or faces for initial review.
- Process one plan at a time.
- Write to a local file geodatabase and confirm adequate free space.
- Disable **Load Results to Map** in the batch tool for large runs.

## A project folder finds no HDF files

**Organize HEC-RAS Project** searches only the selected folder's top level and expects two-digit
names such as `Model.p01.hdf` or `Model.g01.hdf`. Select the HDF directly if it is nested or uses a
different name.

## Layers draw in the wrong location

- Confirm the HEC-RAS project projection and model units.
- Confirm the output feature class CRS.
- Check whether an incorrect override was assigned.
- Compare against known control data rather than relying only on a web basemap.
- Do not use **Define Projection** or another CRS merely to force visual alignment without verifying
  the source coordinates.

## Output names collide with an earlier run

The tools can reuse an existing geodatabase and plan feature dataset. Use a new geodatabase name or
archive the earlier output before rerunning when prior results must be preserved.

## Report a reproducible issue

Open an issue in the [RAS Commander Hydro repository](https://github.com/gpt-cmdr/ras-commander-hydro/issues)
with:

- ArcGIS Pro and Arc Hydro versions.
- HEC-RAS version and model type (1D, 2D, pipe network).
- Tool name and parameter choices.
- Full geoprocessing messages.
- A description of the relevant HDF file type and available groups.
- A small public/reproducible project when licensing and project confidentiality permit.

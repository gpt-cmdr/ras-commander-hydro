# Validation Checklist

Complete these checks before using the extracted layers for engineering interpretation, figures,
spatial analysis, or delivery. The checklist verifies the mechanics of the extraction; it does not
replace project-specific modeling standards or engineering judgment.

## Source control

- [ ] Record the source HEC-RAS project path, geometry number, plan number, and HDF timestamp.
- [ ] Confirm that the selected plan is the intended alternative/scenario.
- [ ] Confirm that computed results are current and that the HEC-RAS compute completed successfully.
- [ ] Keep the original HEC-RAS project unchanged or work from a controlled copy.

## Geoprocessing record

- [ ] Save or capture the ArcGIS Pro geoprocessing messages.
- [ ] Review every warning; do not assume a completed tool produced every selected layer.
- [ ] Record the toolbox/source revision when reproducibility matters.
- [ ] Record the output geodatabase path and whether it existed before the run.

## Spatial reference

- [ ] Confirm the output CRS name, horizontal units, and datum against the HEC-RAS project.
- [ ] Overlay geometry on trusted control data or the project basemap.
- [ ] Check the full model extent for a uniform shift, scale error, or axis/unit mismatch.
- [ ] If **Override CRS** was used, document the authoritative source for that CRS.

## Geometry

- [ ] Compare the number of 1D reaches, cross sections, structures, and 2D flow areas with HEC-RAS.
- [ ] Compare a sample of river/reach names, stations, mesh names, cell IDs, and structure names.
- [ ] Check line direction/orientation where it affects downstream work.
- [ ] Inspect mesh-cell polygons for gaps, overlaps, or malformed cells before polygon analysis.
- [ ] For pipe networks, compare node/conduit counts and key elevations with the source model.

## Results

- [ ] Confirm that result units match the HEC-RAS model unit system.
- [ ] Compare sample maximum WSE and velocity values with RAS Mapper.
- [ ] Compare sample times of maximum with the simulation window.
- [ ] Screen for source NoData/sentinel values before statistics or symbology.
- [ ] Confirm that WSE and velocity layers correspond to the same intended plan.

## Terrain

- [ ] Confirm the `.rasmap` terrain name and VRT path.
- [ ] Check raster extent, resolution, horizontal units, and vertical units.
- [ ] Document that the loaded VRT excludes RAS Mapper vector terrain modifications.
- [ ] If modified terrain is required, compare against the authoritative RAS Mapper terrain workflow.

## Delivery

- [ ] Remove or clearly label temporary/in-memory layers.
- [ ] Use unambiguous plan/scenario names in the geodatabase and map.
- [ ] Include the source/model identification and limitations with figures or shared data.
- [ ] Open the delivered geodatabase in a clean ArcGIS Pro session and confirm that layers load.

Resolve any failed check using the [troubleshooting guide](../troubleshooting.md) or the source
HEC-RAS model before continuing.

# Installation

RAS Commander Arc Hydro Tools runs inside ArcGIS Pro on Windows. It is not an ArcMap toolbox.

## Prerequisites

| Requirement | Notes |
| --- | --- |
| ArcGIS Pro | Use a supported ArcGIS Pro release with a valid license. |
| Arc Hydro Tools | Recommended for the packaged installation. Download the release that matches your ArcGIS Pro version. |
| HEC-RAS project data | The import tools are designed for HEC-RAS 6.x HDF layouts. Results tools require a computed plan HDF. |
| Writable output location | A local file geodatabase is recommended, especially for large meshes. |

## Install with Arc Hydro Tools

This is the recommended path for most users.

1. Download the current Arc Hydro build for your ArcGIS Pro version from the
   [official Esri Arc Hydro download page](https://www.esri.com/en-us/industries/water-resources/arc-hydro/downloads).
2. If another Arc Hydro version is installed, uninstall it first. The Arc Hydro installer does not
   overwrite an existing installation; follow Esri's
   [installation guide](https://www.esri.com/content/dam/esrisites/en-us/media/manuals/downloading-arc-hydro.pdf).
3. Install the downloaded Arc Hydro package with administrator privileges.
4. Restart ArcGIS Pro.
5. In the **Catalog** pane, expand **Toolboxes**, then locate the Arc Hydro toolbox and
   **RAS Commander Tools**.

Arc Hydro packaging can vary by ArcGIS Pro release. If the RAS Commander tools are not present in
the compatible Arc Hydro package available to you, use the source-toolbox method below.

## Add the source toolbox

Use this method for development, evaluation, or access to the latest public repository version.
It does not copy files into the ArcGIS Pro installation directory.

1. Clone or download the repository:

   ```powershell
   git clone https://github.com/gpt-cmdr/ras-commander-hydro.git
   cd ras-commander-hydro
   ```

2. Open ArcGIS Pro and open or create a project with a map.
3. In the **Catalog** pane, right-click **Toolboxes** and select **Add Toolbox**.
4. Browse to `toolboxes/RAS-Commander.pyt` in the repository.
5. Expand **RAS Commander Tools** and open a tool.

Keep the repository structure intact. The Python toolbox loads the supporting `rc_*.py` modules
from `Scripts/archydro` by their relative location.

## Optional system-copy installer

The repository includes `Resources/install_toolbox.ps1` for maintainers who need to copy the
toolbox, scripts, templates, and images into the ArcGIS Pro installation. It requires an elevated
PowerShell session and writes beneath `C:\Program Files\ArcGIS\Pro`. The source-toolbox method is
usually easier to update and remove.

## Verify the installation

Open **Load HEC-RAS 1D Geometry Layers**. A correct installation shows these inputs:

- **Geometry or Plan HDF File**
- **Override CRS (Optional)**
- **Geometry Elements to Load**
- Output feature classes and geodatabase options

If the toolbox has a red error icon, confirm that:

- `RAS-Commander.pyt` still has access to `Scripts/archydro`.
- ArcGIS Pro's Python environment includes `h5py` and `numpy`.
- You installed the Arc Hydro build for your ArcGIS Pro version.

Continue with the [quick start](quick-start.md).

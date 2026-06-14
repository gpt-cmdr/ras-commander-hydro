# Installation

The RAS Commander toolbox is tested in ArcGIS Pro (not the older ArcMap).

## Primary Method: Install with Arc Hydro Tools

The RAS Commander toolbox is included as part of the Arc Hydro Tools distribution. This is the
recommended installation method for most users.

1. Install Arc Hydro Tools with the
   [Arc Hydro installer](https://www.esri.com/en-us/industries/water-resources/arc-hydro/downloads).
2. If you are updating, uninstall Arc Hydro and re-install the latest version (Version 3.4.30
   minimum).
3. The RAS Commander toolbox will be available under:

   ```
   Toolboxes → Arc Hydro Tools → RAS Commander
   ```

## Development Installation

For developers and users who want to extend or customize the tools, or to get the latest
bleeding-edge version between Arc Hydro updates.

1. **Clone the repository**

   ```bash
   git clone https://github.com/gpt-cmdr/ras-commander-hydro.git
   cd ras-commander-hydro
   ```

2. **Option A: Add the toolbox in ArcGIS Pro**

   - Open ArcGIS Pro.
   - In the Catalog pane, right-click on **Toolboxes**.
   - Select **Add Toolbox**.
   - Navigate to `toolboxes/RAS-Commander.pyt`.

3. **Option B: Install for development (requires Admin)**

   ```powershell
   # Run PowerShell as Administrator
   cd Resources
   .\install_toolbox.ps1
   ```

   To uninstall:

   ```powershell
   # Run PowerShell as Administrator
   cd Resources
   .\uninstall_toolbox.ps1
   ```

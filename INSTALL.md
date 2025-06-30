# Installation Instructions for Arc Hydro RAS Commander Toolbox

## Prerequisites

- ArcGIS Pro must be installed
- Administrator privileges are required to install to Program Files

## Installation Methods

### Method 1: Using the Batch File (Easiest)

1. Double-click `install_toolbox_as_admin.bat`
2. Click "Yes" when Windows asks for administrator privileges
3. Follow the prompts in the PowerShell window

### Method 2: Using PowerShell Directly

1. Right-click on `install_toolbox.ps1`
2. Select "Run with PowerShell"
3. If prompted for administrator privileges, click "Yes"

OR

1. Open PowerShell as Administrator
2. Navigate to the repository directory
3. Run: `.\install_toolbox.ps1`

### Method 3: Manual Installation

If the scripts don't work, you can manually copy the files:

1. Copy `Scripts\ras_commander\*` to:
   `C:\Program Files\ArcGIS\Pro\Resources\ArcToolBox\Scripts\ras_commander\`

2. Copy `toolboxes\RAS Commander.pyt` to:
   `C:\Program Files\ArcGIS\Pro\Resources\ArcToolBox\toolboxes\RAS Commander.pyt`

3. (Optional) Copy layer templates, images, and geodatabase templates to their respective directories

## Development Mode

The installation script offers a "development mode" option that creates symbolic links instead of copying files. This is useful for developers because:

- Changes made in the repository are immediately reflected in ArcGIS Pro
- No need to reinstall after making changes
- Easy to switch between development and production versions

To use development mode:
1. Run the installation script
2. When prompted, type 'y' to create symlinks

Note: Creating symbolic links requires administrator privileges on Windows.

## Verifying Installation

After installation:

1. Open ArcGIS Pro
2. Go to the Catalog pane
3. Navigate to Toolboxes → System Toolboxes
4. Look for "RAS Commander"
5. The tools should be available:
   - Load HEC-RAS Terrain
   - Load HEC-RAS 6.x HDF Data

## Troubleshooting

If you encounter issues:

1. **"Script cannot be loaded because running scripts is disabled"**
   - This is a PowerShell execution policy issue
   - The batch file should bypass this, but if not, run:
     `Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser`

2. **"Access denied" errors**
   - Ensure you're running as administrator
   - Check that ArcGIS Pro is not running

3. **"ArcGIS Pro not found"**
   - Verify ArcGIS Pro is installed
   - Check if it's installed in a non-standard location

## Uninstallation

To remove the toolbox:

1. Delete the following:
   - Directory: `C:\Program Files\ArcGIS\Pro\Resources\ArcToolBox\Scripts\ras_commander`
   - File: `C:\Program Files\ArcGIS\Pro\Resources\ArcToolBox\toolboxes\RAS Commander.pyt`

2. (Optional) Remove any layer templates, images, or geodatabase templates that were installed
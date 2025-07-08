# uninstall_toolbox.ps1
#
# Uninstalls the Arc Hydro RAS Commander toolbox from ArcGIS Pro.
# This script removes all toolbox files from the ArcGIS Pro installation directories.
#
# Usage:
#   Right-click on this file and select "Run with PowerShell" 
#   OR
#   Open PowerShell as Administrator and run: .\uninstall_toolbox.ps1
#
# Note: This script requires administrator privileges to remove files from Program Files.

# Check if running as administrator
if (-NOT ([Security.Principal.WindowsPrincipal] [Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole([Security.Principal.WindowsBuiltInRole] "Administrator")) {
    Write-Host "This script requires Administrator privileges." -ForegroundColor Red
    Write-Host "Please run PowerShell as Administrator and try again." -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Press any key to exit..."
    $null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")
    exit 1
}

Write-Host "Arc Hydro RAS Commander Toolbox Uninstaller" -ForegroundColor Cyan
Write-Host ("=" * 50) -ForegroundColor Cyan
Write-Host ""

# Get the script directory (Resources folder) and then get parent (repository root)
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Definition
$repoRoot = Split-Path -Parent $scriptDir

Write-Host "Script location: $scriptDir" -ForegroundColor Gray
Write-Host "Repository root: $repoRoot" -ForegroundColor Gray
Write-Host ""

# Function to find ArcGIS Pro installation
function Find-ArcGISPro {
    $potentialPaths = @(
        "C:\Program Files\ArcGIS\Pro",
        "C:\Program Files (x86)\ArcGIS\Pro",
        "$env:ProgramFiles\ArcGIS\Pro",
        "${env:ProgramFiles(x86)}\ArcGIS\Pro"
    )
    
    # Check each potential path
    foreach ($path in $potentialPaths) {
        if (Test-Path $path) {
            $toolboxPath = Join-Path $path "Resources\ArcToolBox"
            if (Test-Path $toolboxPath) {
                return $path
            }
        }
    }
    
    # Check registry
    try {
        $regPath = "HKLM:\SOFTWARE\ESRI\ArcGISPro"
        if (Test-Path $regPath) {
            $installDir = (Get-ItemProperty -Path $regPath -Name InstallDir -ErrorAction SilentlyContinue).InstallDir
            if ($installDir -and (Test-Path $installDir)) {
                return $installDir
            }
        }
    } catch {
        # Registry check failed, continue
    }
    
    return $null
}

# Find ArcGIS Pro installation
$arcgisProPath = Find-ArcGISPro
if (-not $arcgisProPath) {
    Write-Host "ERROR: Could not find ArcGIS Pro installation." -ForegroundColor Red
    Write-Host "ArcGIS Pro may not be installed or the toolbox may already be uninstalled." -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Press any key to exit..."
    $null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")
    exit 1
}

Write-Host "Found ArcGIS Pro at: $arcgisProPath" -ForegroundColor Green
Write-Host ""

# Define paths to remove
$pathsToRemove = @(
    @{
        Name = "Python Scripts"
        Path = Join-Path $arcgisProPath "Resources\ArcToolBox\Scripts\archydro\rc_*.py"
        Type = "Files"
    },
    @{
        Name = "Python Toolbox (New)"
        Path = Join-Path $arcgisProPath "Resources\ArcToolBox\toolboxes\RAS-Commander.pyt"
        Type = "File"
    },
    @{
        Name = "Python Toolbox (Old)"
        Path = Join-Path $arcgisProPath "Resources\ArcToolBox\toolboxes\RAS Commander.pyt"
        Type = "File"
    },
    @{
        Name = "Python Toolbox XML (New)"
        Path = Join-Path $arcgisProPath "Resources\ArcToolBox\toolboxes\RAS-Commander.pyt.xml"
        Type = "File"
    },
    @{
        Name = "Python Toolbox XML (Old)"
        Path = Join-Path $arcgisProPath "Resources\ArcToolBox\toolboxes\RAS Commander.pyt.xml"
        Type = "File"
    },
    @{
        Name = "Tool-specific XML files"
        Path = Join-Path $arcgisProPath "Resources\ArcToolBox\toolboxes\RAS*.pyt.xml"
        Type = "Files"
    },
    @{
        Name = "Layer Templates"
        Path = Join-Path $arcgisProPath "Resources\ArcToolBox\Templates\Layers\archydro"
        Type = "Directory"
    },
    @{
        Name = "Geodatabase Template"
        Path = Join-Path $arcgisProPath "Resources\ArcToolBox\Data\archydro\Ras2DTemplate.gdb"
        Type = "Directory"
    }
)

# Note: We don't remove Images directory as it might contain files from other tools

Write-Host "This will remove the following Arc Hydro RAS Commander components:" -ForegroundColor Yellow
foreach ($item in $pathsToRemove) {
    if (Test-Path $item.Path) {
        # Check if it's a symlink
        $itemInfo = Get-Item $item.Path -Force -ErrorAction SilentlyContinue
        if ($itemInfo -and $itemInfo.LinkType) {
            Write-Host "  - $($item.Name) (symlink): $($item.Path)" -ForegroundColor Cyan
        } else {
            Write-Host "  - $($item.Name): $($item.Path)" -ForegroundColor White
        }
    }
}

Write-Host ""
$response = Read-Host "Do you want to proceed with uninstallation? (y/N)"

if ($response -ne 'y' -and $response -ne 'Y') {
    Write-Host ""
    Write-Host "Uninstallation cancelled." -ForegroundColor Yellow
    Write-Host "Press any key to exit..."
    $null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")
    exit 0
}

Write-Host ""
Write-Host "Uninstalling Arc Hydro RAS Commander components..." -ForegroundColor Yellow
Write-Host ""

$successCount = 0
$errorCount = 0
$notFoundCount = 0

foreach ($item in $pathsToRemove) {
    Write-Host "Removing: $($item.Name)" -ForegroundColor White
    
    if (-not (Test-Path $item.Path)) {
        Write-Host "  Not found - skipping" -ForegroundColor Gray
        $notFoundCount++
        continue
    }
    
    try {
        # Check if it's a symlink
        $itemInfo = Get-Item $item.Path -Force -ErrorAction SilentlyContinue
        $isSymlink = $itemInfo -and $itemInfo.LinkType
        
        if ($item.Type -eq "Directory") {
            if ($isSymlink) {
                # For symlinked directories, just remove the link
                (Get-Item $item.Path).Delete()
                Write-Host "  Success: Removed symlink" -ForegroundColor Green
            } else {
                # For regular directories, remove recursively
                Remove-Item -Path $item.Path -Recurse -Force
                Write-Host "  Success: Removed directory and all contents" -ForegroundColor Green
            }
        } elseif ($item.Type -eq "Files") {
            # For file patterns (like rc_*.py)
            $parentPath = Split-Path $item.Path -Parent
            $pattern = Split-Path $item.Path -Leaf
            Get-ChildItem -Path $parentPath -Filter $pattern -ErrorAction SilentlyContinue | ForEach-Object {
                Remove-Item -Path $_.FullName -Force
            }
            Write-Host "  Success: Removed matching files" -ForegroundColor Green
        } else {
            # For single files (including symlinked files)
            Remove-Item -Path $item.Path -Force
            Write-Host "  Success: Removed file" -ForegroundColor Green
        }
        
        $successCount++
        
    } catch {
        Write-Host "  ERROR: Failed to remove - $_" -ForegroundColor Red
        $errorCount++
    }
}

Write-Host ""
Write-Host ("=" * 50) -ForegroundColor Cyan
Write-Host "Uninstallation Summary:" -ForegroundColor Cyan
Write-Host "  Removed: $successCount components" -ForegroundColor Green
Write-Host "  Failed: $errorCount components" -ForegroundColor $(if ($errorCount -gt 0) { "Red" } else { "Gray" })
Write-Host "  Not found: $notFoundCount components" -ForegroundColor Gray
Write-Host ""

if ($errorCount -gt 0) {
    Write-Host "Some components failed to uninstall. Please check the errors above." -ForegroundColor Yellow
    Write-Host "You may need to manually remove these items." -ForegroundColor Yellow
} elseif ($successCount -eq 0 -and $notFoundCount -gt 0) {
    Write-Host "No components were found to uninstall." -ForegroundColor Yellow
    Write-Host "The toolbox may have already been uninstalled." -ForegroundColor Yellow
} else {
    Write-Host "Arc Hydro RAS Commander has been successfully uninstalled!" -ForegroundColor Green
}

# Clean up empty parent directories if they exist
Write-Host ""
Write-Host "Cleaning up empty directories..." -ForegroundColor Yellow

$parentDirsToCheck = @(
    Join-Path $arcgisProPath "Resources\ArcToolBox\Scripts\ras_commander"
    Join-Path $arcgisProPath "Resources\ArcToolBox\Templates\Layers\archydro"
    Join-Path $arcgisProPath "Resources\ArcToolBox\Data\archydro"
)

foreach ($dir in $parentDirsToCheck) {
    $parent = Split-Path $dir -Parent
    if ((Test-Path $parent) -and (Get-ChildItem $parent -Force | Measure-Object).Count -eq 0) {
        try {
            Remove-Item $parent -Force
            Write-Host "  Removed empty directory: $parent" -ForegroundColor Green
        } catch {
            # Ignore errors for cleanup
        }
    }
}

Write-Host ""
Write-Host "Press any key to exit..."
$null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")
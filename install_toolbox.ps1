# install_toolbox.ps1
#
# Installs the Arc Hydro RAS Commander toolbox for local development.
# This script copies all toolbox files from the repository to the ArcGIS Pro installation directories.
#
# Usage:
#   Right-click on this file and select "Run with PowerShell" 
#   OR
#   Open PowerShell as Administrator and run: .\install_toolbox.ps1
#
# Note: This script requires administrator privileges to write to Program Files.

# Check if running as administrator
if (-NOT ([Security.Principal.WindowsPrincipal] [Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole([Security.Principal.WindowsBuiltInRole] "Administrator")) {
    Write-Host "This script requires Administrator privileges." -ForegroundColor Red
    Write-Host "Please run PowerShell as Administrator and try again." -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Press any key to exit..."
    $null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")
    exit 1
}

Write-Host "Arc Hydro RAS Commander Toolbox Installer" -ForegroundColor Cyan
Write-Host ("=" * 50) -ForegroundColor Cyan
Write-Host ""

# Get the script directory (repository root)
$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Definition

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
    Write-Host "Please ensure ArcGIS Pro is installed." -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Press any key to exit..."
    $null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")
    exit 1
}

Write-Host "Found ArcGIS Pro at: $arcgisProPath" -ForegroundColor Green
Write-Host ""

# Define source and destination mappings
$mappings = @(
    @{
        Name = "Python Scripts"
        Source = Join-Path $repoRoot "Scripts\archydro"
        Destination = Join-Path $arcgisProPath "Resources\ArcToolBox\Scripts\archydro"
        Type = "Directory"
        Required = $true
        Filter = "rc_*.py"
    },
    @{
        Name = "Python Toolbox"
        Source = Join-Path $repoRoot "toolboxes\RAS Commander.pyt"
        Destination = Join-Path $arcgisProPath "Resources\ArcToolBox\toolboxes\RAS Commander.pyt"
        Type = "File"
        Required = $true
    },
    @{
        Name = "Layer Templates"
        Source = Join-Path $repoRoot "Templates\Layers\archydro"
        Destination = Join-Path $arcgisProPath "Resources\ArcToolBox\Templates\Layers\archydro"
        Type = "Directory"
        Required = $false
    },
    @{
        Name = "Images"
        Source = Join-Path $repoRoot "Images"
        Destination = Join-Path $arcgisProPath "Resources\ArcToolBox\Images"
        Type = "Directory"
        Required = $false
    },
    @{
        Name = "Geodatabase Template"
        Source = Join-Path $repoRoot "Data\archydro\Ras2DTemplate.gdb"
        Destination = Join-Path $arcgisProPath "Resources\ArcToolBox\Data\archydro\Ras2DTemplate.gdb"
        Type = "Directory"
        Required = $false
    }
)

Write-Host "Installing Arc Hydro RAS Commander components..." -ForegroundColor Yellow
Write-Host ""

$successCount = 0
$errorCount = 0
$skippedCount = 0

foreach ($mapping in $mappings) {
    Write-Host "Installing: $($mapping.Name)" -ForegroundColor White
    
    if (-not (Test-Path $mapping.Source)) {
        if ($mapping.Required) {
            Write-Host "  ERROR: Source not found: $($mapping.Source)" -ForegroundColor Red
            $errorCount++
        } else {
            Write-Host "  Skipping (optional component not found)" -ForegroundColor Gray
            $skippedCount++
        }
        continue
    }
    
    try {
        # Create parent directory if it doesn't exist
        $parentDir = Split-Path -Parent $mapping.Destination
        if (-not (Test-Path $parentDir)) {
            New-Item -ItemType Directory -Path $parentDir -Force | Out-Null
        }
        
        # Remove existing destination if it exists
        if (Test-Path $mapping.Destination) {
            Remove-Item -Path $mapping.Destination -Recurse -Force
        }
        
        # Copy the content
        if ($mapping.Type -eq "Directory") {
            if ($mapping.Filter) {
                # Create destination directory
                New-Item -ItemType Directory -Path $mapping.Destination -Force | Out-Null
                # Copy only filtered files
                Get-ChildItem -Path $mapping.Source -Filter $mapping.Filter | ForEach-Object {
                    Copy-Item -Path $_.FullName -Destination $mapping.Destination -Force
                }
            } else {
                Copy-Item -Path $mapping.Source -Destination $mapping.Destination -Recurse -Force
            }
        } else {
            # For files, ensure the destination directory exists
            $destDir = Split-Path -Parent $mapping.Destination
            if (-not (Test-Path $destDir)) {
                New-Item -ItemType Directory -Path $destDir -Force | Out-Null
            }
            Copy-Item -Path $mapping.Source -Destination $mapping.Destination -Force
        }
        
        Write-Host "  Success: Copied to $($mapping.Destination)" -ForegroundColor Green
        $successCount++
        
    } catch {
        Write-Host "  ERROR: Failed to copy - $_" -ForegroundColor Red
        $errorCount++
    }
}

Write-Host ""
Write-Host ("=" * 50) -ForegroundColor Cyan
Write-Host "Installation Summary:" -ForegroundColor Cyan
Write-Host "  Successful: $successCount components" -ForegroundColor Green
Write-Host "  Failed: $errorCount components" -ForegroundColor $(if ($errorCount -gt 0) { "Red" } else { "Gray" })
Write-Host "  Skipped: $skippedCount optional components" -ForegroundColor Gray
Write-Host ""

if ($errorCount -gt 0) {
    Write-Host "Some components failed to install. Please check the errors above." -ForegroundColor Yellow
} else {
    Write-Host "All components installed successfully!" -ForegroundColor Green
    Write-Host ""
    Write-Host "To use the toolbox in ArcGIS Pro:" -ForegroundColor Yellow
    Write-Host "  1. Open ArcGIS Pro"
    Write-Host "  2. Go to the Catalog pane"
    Write-Host "  3. Navigate to Toolboxes > Arc Hydro RAS Commander"
    Write-Host "  4. The tools will be available there"
}

Write-Host ""
Write-Host ("-" * 50) -ForegroundColor Gray

# Ask about development mode (symlinks)
Write-Host ""
$response = Read-Host "Would you like to create development symlinks instead of copying? (y/N)"

if ($response -eq 'y' -or $response -eq 'Y') {
    Write-Host ""
    Write-Host "Creating symlinks for development mode..." -ForegroundColor Yellow
    Write-Host "(This allows changes in the repo to be reflected immediately)" -ForegroundColor Gray
    Write-Host ""
    
    $symlinkSuccess = 0
    $symlinkError = 0
    
    foreach ($mapping in $mappings) {
        if (-not (Test-Path $mapping.Source)) {
            continue
        }
        
        Write-Host "Creating symlink: $($mapping.Name)" -ForegroundColor White
        
        try {
            # Remove existing destination if it exists
            if (Test-Path $mapping.Destination) {
                Remove-Item -Path $mapping.Destination -Recurse -Force
            }
            
            # Create parent directory if needed
            $parentDir = Split-Path -Parent $mapping.Destination
            if (-not (Test-Path $parentDir)) {
                New-Item -ItemType Directory -Path $parentDir -Force | Out-Null
            }
            
            # Create symlink
            New-Item -ItemType SymbolicLink -Path $mapping.Destination -Target $mapping.Source | Out-Null
            
            Write-Host "  Success: Symlinked to $($mapping.Destination)" -ForegroundColor Green
            $symlinkSuccess++
            
        } catch {
            Write-Host "  ERROR: Failed to create symlink - $_" -ForegroundColor Red
            $symlinkError++
        }
    }
    
    Write-Host ""
    Write-Host "Symlink Summary:" -ForegroundColor Cyan
    Write-Host "  Successful: $symlinkSuccess symlinks" -ForegroundColor Green
    Write-Host "  Failed: $symlinkError symlinks" -ForegroundColor $(if ($symlinkError -gt 0) { "Red" } else { "Gray" })
}

Write-Host ""
Write-Host "Press any key to exit..."
$null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")
@echo off
REM update_metadata.bat - Updates all RAS Commander tool metadata XML files

echo RAS Commander Metadata Updater
echo ==============================
echo.

REM Set path to ArcGIS Pro Python
set ARCGIS_PYTHON="C:\Program Files\ArcGIS\Pro\bin\Python\envs\arcgispro-py3\python.exe"

REM Check if ArcGIS Pro Python exists
if not exist %ARCGIS_PYTHON% (
    echo ERROR: ArcGIS Pro Python not found at %ARCGIS_PYTHON%
    echo Please update the path in this batch file to match your installation
    pause
    exit /b 1
)

echo Using ArcGIS Pro Python: %ARCGIS_PYTHON%
echo.

echo Updating tool metadata XML files...
echo.

REM Run the metadata update script with ArcGIS Pro Python
%ARCGIS_PYTHON% update_tool_metadata.py

if %errorlevel% equ 0 (
    echo.
    echo Metadata update completed successfully!
    echo.
    echo Updated files:
    echo - toolboxes\RAS Commander.LoadHECRAS1DGeometry.pyt.xml
    echo - toolboxes\RAS Commander.LoadHECRAS2DGeometry.pyt.xml
    echo - toolboxes\RAS Commander.LoadHECRAS2DResults.pyt.xml
    echo - toolboxes\RAS Commander.LoadRASTerrain.pyt.xml
    echo - toolboxes\RAS Commander.OrganizeRASProject.pyt.xml
    echo - toolboxes\RAS Commander.pyt.xml
) else (
    echo.
    echo ERROR: Metadata update failed!
    echo Check the error messages above for details.
)

echo.
pause
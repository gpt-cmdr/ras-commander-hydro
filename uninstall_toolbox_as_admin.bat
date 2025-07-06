@echo off
:: uninstall_toolbox_as_admin.bat
:: Wrapper to run the uninstall PowerShell script with administrator privileges

echo Arc Hydro RAS Commander Toolbox Uninstaller
echo ============================================
echo.
echo This will uninstall the RAS Commander toolbox from ArcGIS Pro.
echo Administrator privileges are required.
echo.

:: Run PowerShell script with admin privileges and bypass execution policy
powershell.exe -ExecutionPolicy Bypass -Command "Start-Process PowerShell -ArgumentList '-ExecutionPolicy Bypass -File ""%~dp0uninstall_toolbox.ps1""' -Verb RunAs"
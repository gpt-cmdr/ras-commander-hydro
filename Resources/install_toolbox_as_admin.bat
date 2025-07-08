@echo off
:: install_toolbox_as_admin.bat
:: 
:: Convenience batch file to run the PowerShell installation script as administrator
:: Just double-click this file to install the Arc Hydro RAS Commander toolbox

echo Arc Hydro RAS Commander Toolbox Installer
echo =========================================
echo.
echo This will install the toolbox to your ArcGIS Pro installation.
echo Administrator privileges are required.
echo.
pause

:: Run PowerShell script as administrator
powershell -Command "Start-Process powershell -ArgumentList '-NoProfile -ExecutionPolicy Bypass -File ""%~dp0install_toolbox.ps1""' -Verb RunAs"
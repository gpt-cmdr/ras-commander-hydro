@echo off
REM ============================================================================
REM RAS Commander Arc Hydro Tools - Test Runner
REM
REM Usage:
REM   run_tests.bat              - Run tier1 tests (default)
REM   run_tests.bat tier1        - Run tier1 tests only
REM   run_tests.bat tier2        - Run tier2 (production model) tests only
REM   run_tests.bat all          - Run all tests
REM   run_tests.bat unicode      - Run unicode-specific tests
REM   run_tests.bat unit         - Run unit tests only
REM   run_tests.bat integration  - Run integration tests only
REM
REM Additional pytest args can be appended after the mode:
REM   run_tests.bat tier1 -x --pdb
REM ============================================================================

setlocal enabledelayedexpansion

REM Find propy.bat - ArcGIS Pro's Python launcher
set "PROPY="

REM Check common installation paths
for %%P in (
    "%ProgramFiles%\ArcGIS\Pro\bin\Python\Scripts\propy.bat"
    "%ProgramFiles(x86)%\ArcGIS\Pro\bin\Python\Scripts\propy.bat"
    "C:\Program Files\ArcGIS\Pro\bin\Python\Scripts\propy.bat"
    "%LOCALAPPDATA%\Programs\ArcGIS\Pro\bin\Python\Scripts\propy.bat"
) do (
    if exist %%P (
        set "PROPY=%%~P"
        goto :found_propy
    )
)

REM Try to find via registry
for /f "tokens=2*" %%A in ('reg query "HKLM\SOFTWARE\ESRI\ArcGISPro" /v InstallDir 2^>nul') do (
    set "ARCGIS_DIR=%%B"
    if exist "!ARCGIS_DIR!bin\Python\Scripts\propy.bat" (
        set "PROPY=!ARCGIS_DIR!bin\Python\Scripts\propy.bat"
        goto :found_propy
    )
)

echo ERROR: Could not find ArcGIS Pro's propy.bat
echo Please ensure ArcGIS Pro is installed, or set the PROPY environment variable.
echo Example: set PROPY=C:\Program Files\ArcGIS\Pro\bin\Python\Scripts\propy.bat
exit /b 1

:found_propy
echo Using Python: %PROPY%

REM Parse the test mode argument
set "MODE=%~1"
set "EXTRA_ARGS="

REM Collect extra args (everything after the first argument)
if "%MODE%"=="" set "MODE=tier1"
shift
:collect_args
if "%~1"=="" goto :run_tests
set "EXTRA_ARGS=%EXTRA_ARGS% %~1"
shift
goto :collect_args

:run_tests
REM Map mode to pytest marker expression
if /i "%MODE%"=="tier1" (
    echo Running Tier 1 tests...
    "%PROPY%" -m pytest tests/ -m "tier1 and not tier2" %EXTRA_ARGS%
) else if /i "%MODE%"=="tier2" (
    echo Running Tier 2 ^(production model^) tests...
    "%PROPY%" -m pytest tests/tier2/ -m "tier2" %EXTRA_ARGS%
) else if /i "%MODE%"=="all" (
    echo Running all tests...
    "%PROPY%" -m pytest tests/ %EXTRA_ARGS%
) else if /i "%MODE%"=="unicode" (
    echo Running unicode tests...
    "%PROPY%" -m pytest tests/ -m "unicode" %EXTRA_ARGS%
) else if /i "%MODE%"=="unit" (
    echo Running unit tests...
    "%PROPY%" -m pytest tests/unit/ %EXTRA_ARGS%
) else if /i "%MODE%"=="integration" (
    echo Running integration tests...
    "%PROPY%" -m pytest tests/integration/ %EXTRA_ARGS%
) else (
    echo Unknown mode: %MODE%
    echo Valid modes: tier1, tier2, all, unicode, unit, integration
    exit /b 1
)

exit /b %ERRORLEVEL%

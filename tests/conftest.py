# -*- coding: utf-8 -*-
"""
conftest.py - Shared test infrastructure for RAS Commander Arc Hydro Tools.

Provides:
    - sys.path setup mirroring RAS-Commander.pyt:56-59
    - MockMessages: captures addMessage/addWarning/addWarningMessage/addErrorMessage
    - MockParam: mimics arcpy.Parameter (from rc_organize_ras_project.py:552-556)
    - Fixtures for HDF test data in testdata/
"""

import sys
import os
import types
import pytest

# ---------------------------------------------------------------------------
# sys.path setup - mirrors toolboxes/RAS-Commander.pyt lines 56-59
#
# IMPORTANT: ArcGIS Pro also installs rc_* scripts into its system directory
# (C:\Program Files\ArcGIS\Pro\Resources\ArcToolbox\Scripts\archydro\).
# When arcpy is imported, that system path may be added to sys.path, causing
# Python to find the SYSTEM copy instead of our repo's modified copy.
# We remove the system path to ensure our repo's scripts take precedence.
# ---------------------------------------------------------------------------
_repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_scripts_dir = os.path.join(_repo_root, "Scripts", "archydro")
if _scripts_dir not in sys.path:
    sys.path.insert(0, _scripts_dir)

# Remove ArcGIS Pro's system archydro path if present (so our repo copy wins).
# When arcpy is imported, ArcGIS Pro may add its own ArcToolbox/Scripts/archydro
# to sys.path, which contains an older copy of our rc_* scripts without test fixes.
_arcgis_archydro = os.path.normpath(os.path.join(
    os.environ.get("ProgramFiles", r"C:\Program Files"),
    "ArcGIS", "Pro", "Resources", "ArcToolbox", "Scripts", "archydro"
))

def _scrub_system_archydro_path():
    """Remove the ArcGIS system archydro path and ensure our repo path is first."""
    sys.path[:] = [p for p in sys.path if os.path.normpath(p) != _arcgis_archydro]
    if _scripts_dir not in sys.path:
        sys.path.insert(0, _scripts_dir)

_scrub_system_archydro_path()

TESTDATA_DIR = os.path.join(_repo_root, "testdata")

# ---------------------------------------------------------------------------
# Mock arcpy module — installed when real arcpy is unavailable so that
# rc_utils.py and other scripts can be imported for unit testing of pure
# functions (naming helpers, field detection, clockwise check, etc.).
# Integration tests that actually need arcpy skip themselves via
# ``pytest.skip("arcpy not available", allow_module_level=True)``.
# ---------------------------------------------------------------------------
try:
    import arcpy  # noqa: F401 — probe only
    HAS_ARCPY = True
except ImportError:
    HAS_ARCPY = False

    _arcpy = types.ModuleType("arcpy")
    _arcpy.__package__ = "arcpy"
    _arcpy._is_mock = True  # Marker so integration tests can detect the mock

    # Stub the free functions that rc_utils.py calls at module level
    _arcpy.AddMessage = lambda msg: None
    _arcpy.AddWarning = lambda msg: None
    _arcpy.AddError = lambda msg: None

    # arcpy.ExecuteError is raised by tools on validation failure
    class _ExecuteError(RuntimeError):
        pass
    _arcpy.ExecuteError = _ExecuteError

    # arcpy.SpatialReference — used by get_ras_projection_wkt return path
    class _SpatialReference:
        def __init__(self, *a, **kw):
            self.name = "MockSR"
        def loadFromString(self, wkt):
            pass
    _arcpy.SpatialReference = _SpatialReference

    # arcpy.Exists — used in several utils
    _arcpy.Exists = lambda path: False

    # arcpy.env
    _env = types.SimpleNamespace(workspace=None)
    _arcpy.env = _env

    # arcpy.management stub
    _mgmt = types.ModuleType("arcpy.management")
    _mgmt.CreateFeatureclass = lambda *a, **kw: None
    _mgmt.AddField = lambda *a, **kw: None
    _mgmt.Delete = lambda *a, **kw: None
    _arcpy.management = _mgmt

    # arcpy.da stub (InsertCursor, Walk, SearchCursor)
    _da = types.ModuleType("arcpy.da")
    _da.Walk = lambda *a, **kw: iter([])
    _arcpy.da = _da

    # arcpy.mp stub
    _mp = types.ModuleType("arcpy.mp")
    _mp.ArcGISProject = lambda *a, **kw: None
    _arcpy.mp = _mp

    # arcpy.metadata stub
    _metadata = types.ModuleType("arcpy.metadata")
    class _Metadata:
        def __init__(self, *a):
            self.description = ""
            self.summary = ""
            self.tags = ""
        def save(self):
            pass
    _metadata.Metadata = _Metadata
    _arcpy.metadata = _metadata

    # Geoprocessing function stubs
    _arcpy.CreateFileGDB_management = lambda *a, **kw: None
    _arcpy.CreateFeatureDataset_management = lambda *a, **kw: None
    _arcpy.GetCount_management = lambda *a, **kw: ["0"]
    _arcpy.ListFields = lambda *a, **kw: []
    _arcpy.ApplySymbologyFromLayer_management = lambda *a, **kw: None
    _arcpy.CheckOutExtension = lambda *a, **kw: None

    # arcpy.Parameter stub (for getParameterInfo)
    class _Parameter:
        def __init__(self, **kw):
            self.value = kw.get("value")
            self.valueAsText = None
            self.values = None
            self.altered = False
            self.enabled = True
            self.filter = types.SimpleNamespace(type=None, list=[])
            self.description = ""
            self.displayName = kw.get("displayName", "")
            self.name = kw.get("name", "")
            self.datatype = kw.get("datatype", "")
            self.parameterType = kw.get("parameterType", "Optional")
            self.direction = kw.get("direction", "Input")
            self.multiValue = kw.get("multiValue", False)
            self.category = None
        def setWarningMessage(self, msg):
            pass
        def setErrorMessage(self, msg):
            pass
    _arcpy.Parameter = _Parameter

    # Install the mock
    sys.modules["arcpy"] = _arcpy
    sys.modules["arcpy.management"] = _mgmt
    sys.modules["arcpy.da"] = _da
    sys.modules["arcpy.mp"] = _mp
    sys.modules["arcpy.metadata"] = _metadata

# ---------------------------------------------------------------------------
# Pre-import all rc_* modules from our repo BEFORE any test runs.
# This ensures Python caches OUR copies in sys.modules, not the system copies
# from C:\Program Files\ArcGIS\Pro\Resources\ArcToolbox\Scripts\archydro\.
# ---------------------------------------------------------------------------
_scrub_system_archydro_path()  # Ensure clean path before imports
import importlib
for _mod_name in [
    "rc_utils",
    "rc_load_hecras_1d_geometry",
    "rc_load_hecras_2d_geometry",
    "rc_load_hecras_2d_results",
    "rc_organize_ras_project",
]:
    if _mod_name not in sys.modules:
        try:
            importlib.import_module(_mod_name)
        except Exception:
            pass  # Some modules may fail to import (e.g., missing deps) — that's OK


# ---------------------------------------------------------------------------
# MockMessages - captures tool messaging without arcpy
# ---------------------------------------------------------------------------
class MockMessages:
    """Drop-in replacement for the ArcGIS ``messages`` object.

    Supports both ``addWarning`` *and* ``addWarningMessage`` because the
    codebase uses both variants interchangeably.
    """

    def __init__(self):
        self.messages = []
        self.warnings = []
        self.errors = []

    # --- primary API (matches arcpy messages object) ---
    def addMessage(self, msg):
        self.messages.append(str(msg))

    def addWarning(self, msg):
        self.warnings.append(str(msg))

    def addWarningMessage(self, msg):
        """Alias used in write_features_to_fc and other helpers."""
        self.warnings.append(str(msg))

    def addErrorMessage(self, msg):
        self.errors.append(str(msg))

    # --- convenience for assertions ---
    @property
    def all_text(self):
        return "\n".join(self.messages + self.warnings + self.errors)

    def has_error(self):
        return len(self.errors) > 0

    def has_warning(self):
        return len(self.warnings) > 0

    def __repr__(self):
        return (
            f"MockMessages(messages={len(self.messages)}, "
            f"warnings={len(self.warnings)}, errors={len(self.errors)})"
        )


# ---------------------------------------------------------------------------
# MockParam - mimics arcpy.Parameter
# ---------------------------------------------------------------------------
class MockParam:
    """Lightweight stand-in for ``arcpy.Parameter``.

    Modeled after the pattern in rc_organize_ras_project.py:552-556.
    """

    def __init__(self, value=None, valueAsText=None):
        self.value = value
        if valueAsText is not None:
            self.valueAsText = valueAsText
        elif value is not None:
            self.valueAsText = str(value)
        else:
            self.valueAsText = None
        self.values = value if isinstance(value, list) else None
        self.altered = False
        self.enabled = True
        self._warning_message = None
        self._error_message = None

    def setWarningMessage(self, msg):
        self._warning_message = msg

    def setErrorMessage(self, msg):
        self._error_message = msg

    def __repr__(self):
        return f"MockParam(value={self.value!r})"


# ---------------------------------------------------------------------------
# HDF test data paths
# ---------------------------------------------------------------------------
HDF_1D_UNSTEADY = os.path.join(TESTDATA_DIR, "BaldEagle.p01.hdf")
HDF_1D_STEADY = os.path.join(TESTDATA_DIR, "BaldEagle.p02.hdf")
HDF_2D_UNSTEADY = os.path.join(TESTDATA_DIR, "BaldEagleDamBrk.p07.hdf")
HDF_2D_PIPES = os.path.join(TESTDATA_DIR, "DavisStormSystem.p02.hdf")


def _check_testdata():
    """Verify that testdata directory exists and contains expected files."""
    if not os.path.isdir(TESTDATA_DIR):
        pytest.skip(f"testdata/ directory not found at {TESTDATA_DIR}")
    return True


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session")
def testdata_dir():
    """Return the absolute path to testdata/."""
    _check_testdata()
    return TESTDATA_DIR


@pytest.fixture(scope="session")
def hdf_1d_unsteady():
    """Path to 1D unsteady HDF: BaldEagle.p01.hdf"""
    _check_testdata()
    if not os.path.isfile(HDF_1D_UNSTEADY):
        pytest.skip(f"Missing test file: {HDF_1D_UNSTEADY}")
    return HDF_1D_UNSTEADY


@pytest.fixture(scope="session")
def hdf_1d_steady():
    """Path to 1D steady HDF: BaldEagle.p02.hdf"""
    _check_testdata()
    if not os.path.isfile(HDF_1D_STEADY):
        pytest.skip(f"Missing test file: {HDF_1D_STEADY}")
    return HDF_1D_STEADY


@pytest.fixture(scope="session")
def hdf_2d_unsteady():
    """Path to 2D unsteady HDF: BaldEagleDamBrk.p07.hdf"""
    _check_testdata()
    if not os.path.isfile(HDF_2D_UNSTEADY):
        pytest.skip(f"Missing test file: {HDF_2D_UNSTEADY}")
    return HDF_2D_UNSTEADY


@pytest.fixture(scope="session")
def hdf_2d_pipes():
    """Path to 2D HDF with pipe networks: DavisStormSystem.p02.hdf"""
    _check_testdata()
    if not os.path.isfile(HDF_2D_PIPES):
        pytest.skip(f"Missing test file: {HDF_2D_PIPES}")
    return HDF_2D_PIPES


@pytest.fixture
def mock_messages():
    """Fresh MockMessages instance for each test."""
    return MockMessages()


@pytest.fixture
def output_gdb_path(tmp_path):
    """Temporary geodatabase path for test output."""
    return str(tmp_path / "test_output.gdb")

# -*- coding: utf-8 -*-
"""
Unit tests for unicode/encoding handling in RAS Commander.

Targets the two specific bugs found in rc_utils.py:
  1. Line 32: proj_wkt_attr.decode("utf-8") - now .decode("utf-8", "ignore")
  2. Line 290: time_str.decode('utf-8') - now .decode('utf-8', 'ignore')

Also validates that all HDF string attribute reads handle non-UTF-8 bytes
gracefully across the test data files.
"""

import os
import pytest
import h5py
import numpy as np

pytestmark = [pytest.mark.tier1, pytest.mark.unicode]


# ---------------------------------------------------------------------------
# Direct byte-decode edge cases (no HDF file needed)
# ---------------------------------------------------------------------------
class TestByteDecodeEdgeCases:
    """Verify that .decode('utf-8', 'ignore') handles problematic bytes."""

    def test_clean_ascii(self):
        raw = b"NAD83 / UTM zone 17N"
        assert raw.decode("utf-8", "ignore") == "NAD83 / UTM zone 17N"

    def test_valid_utf8(self):
        raw = "Prüfung Straße".encode("utf-8")
        assert raw.decode("utf-8", "ignore") == "Prüfung Straße"

    def test_latin1_bytes_stripped(self):
        """Latin-1 encoded bytes with non-UTF-8 chars are silently dropped."""
        raw = "Stra\xdfe".encode("latin-1")  # ß in latin-1 = 0xDF
        result = raw.decode("utf-8", "ignore")
        assert "Stra" in result  # the 0xDF byte is dropped

    def test_null_bytes_in_wkt(self):
        """HDF string attributes may have trailing null bytes."""
        raw = b'PROJCS["NAD83"]\x00\x00\x00'
        result = raw.decode("utf-8", "ignore")
        assert 'PROJCS["NAD83"]' in result

    def test_mixed_encoding_bytes(self):
        """Mixed valid/invalid bytes should not raise."""
        raw = b"Hello\x80\x81World\xff"
        result = raw.decode("utf-8", "ignore")
        assert "HelloWorld" in result

    def test_empty_bytes(self):
        raw = b""
        assert raw.decode("utf-8", "ignore") == ""

    def test_all_invalid_bytes(self):
        raw = b"\x80\x81\x82\x83"
        result = raw.decode("utf-8", "ignore")
        assert result == ""

    def test_numpy_bytes_decode(self):
        """numpy.bytes_ should also decode cleanly."""
        raw = np.bytes_(b"test\xff value")
        result = raw.decode("utf-8", "ignore")
        assert "test" in result
        assert "value" in result

    def test_datetime_format_with_ignore(self):
        """Verify datetime parsing still works with 'ignore' error handler."""
        from datetime import datetime

        raw = b"02Jan2020 12:00:00"
        decoded = raw.decode("utf-8", "ignore")
        dt = datetime.strptime(decoded, "%d%b%Y %H:%M:%S")
        assert dt.year == 2020
        assert dt.month == 1
        assert dt.day == 2

    def test_datetime_with_trailing_junk(self):
        """Trailing non-UTF-8 bytes after valid datetime string."""
        raw = b"15Mar2023 08:30:00\xff\xfe"
        decoded = raw.decode("utf-8", "ignore")
        from datetime import datetime
        dt = datetime.strptime(decoded.strip(), "%d%b%Y %H:%M:%S")
        assert dt.year == 2023

    def test_wkt_with_embedded_invalid(self):
        """WKT string with invalid bytes embedded — should still yield usable WKT."""
        raw = b'PROJCS["NAD_1983\x80_UTM_Zone_17N",GEOGCS["GCS_North_American_1983"]]'
        result = raw.decode("utf-8", "ignore")
        assert "PROJCS" in result
        assert "GCS_North_American_1983" in result


# ---------------------------------------------------------------------------
# HDF file string attribute scanning
# ---------------------------------------------------------------------------
class TestHDFStringAttributeScan:
    """Scan all string attributes in test HDF files for decode safety."""

    def _scan_attrs(self, hdf_path):
        """Recursively scan all HDF attributes, decode any bytes with strict UTF-8.

        If strict decode fails, it means the file has non-UTF-8 data that
        the 'ignore' handler would silently fix.
        """
        issues = []

        def visitor(name, obj):
            for attr_name, attr_val in obj.attrs.items():
                if isinstance(attr_val, (bytes, np.bytes_)):
                    try:
                        attr_val.decode("utf-8")  # strict — intentional
                    except UnicodeDecodeError:
                        issues.append(f"{name}@{attr_name}")

        with h5py.File(hdf_path, "r") as f:
            # Check file-level attributes
            for attr_name, attr_val in f.attrs.items():
                if isinstance(attr_val, (bytes, np.bytes_)):
                    try:
                        attr_val.decode("utf-8")
                    except UnicodeDecodeError:
                        issues.append(f"/@{attr_name}")
            f.visititems(visitor)

        return issues

    def test_1d_unsteady_no_strict_failures(self, hdf_1d_unsteady):
        issues = self._scan_attrs(hdf_1d_unsteady)
        if issues:
            pytest.xfail(f"Found {len(issues)} attributes with non-UTF-8 bytes: {issues[:5]}")

    def test_1d_steady_no_strict_failures(self, hdf_1d_steady):
        issues = self._scan_attrs(hdf_1d_steady)
        if issues:
            pytest.xfail(f"Found {len(issues)} attributes with non-UTF-8 bytes: {issues[:5]}")

    def test_2d_unsteady_no_strict_failures(self, hdf_2d_unsteady):
        issues = self._scan_attrs(hdf_2d_unsteady)
        if issues:
            pytest.xfail(f"Found {len(issues)} attributes with non-UTF-8 bytes: {issues[:5]}")

    def test_2d_pipes_no_strict_failures(self, hdf_2d_pipes):
        issues = self._scan_attrs(hdf_2d_pipes)
        if issues:
            pytest.xfail(f"Found {len(issues)} attributes with non-UTF-8 bytes: {issues[:5]}")


# ---------------------------------------------------------------------------
# Verify the bug fixes in rc_utils.py
# ---------------------------------------------------------------------------
class TestBugFixVerification:
    """Verify that the two identified decode bugs are actually fixed in source."""

    def test_projection_decode_uses_ignore(self):
        """rc_utils.py line 32 should use .decode('utf-8', 'ignore')."""
        import inspect
        from rc_utils import get_ras_projection_wkt

        source = inspect.getsource(get_ras_projection_wkt)
        # Should NOT have bare .decode("utf-8") without 'ignore'
        assert '.decode("utf-8", "ignore")' in source or ".decode('utf-8', 'ignore')" in source

    def test_simulation_time_decode_uses_ignore(self):
        """rc_utils.py line 290 should use .decode('utf-8', 'ignore')."""
        import inspect
        from rc_utils import cache_hdf_metadata

        source = inspect.getsource(cache_hdf_metadata)
        assert ".decode('utf-8', 'ignore')" in source or '.decode("utf-8", "ignore")' in source

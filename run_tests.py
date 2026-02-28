#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
RAS Commander Arc Hydro Tools - Python Test Entry Point

Usage:
    python run_tests.py                    # Run tier1 tests (default)
    python run_tests.py tier1              # Run tier1 tests
    python run_tests.py tier2              # Run tier2 (production) tests
    python run_tests.py all                # Run all tests
    python run_tests.py unicode            # Run unicode tests
    python run_tests.py unit               # Run unit tests only
    python run_tests.py -- -x --pdb        # Pass extra pytest args after --
"""

import sys
import subprocess


def main():
    args = sys.argv[1:]

    # Separate mode from extra pytest args
    mode = "tier1"
    extra_args = []

    if args:
        if args[0] == "--":
            extra_args = args[1:]
        elif args[0] in ("tier1", "tier2", "all", "unicode", "unit", "integration"):
            mode = args[0]
            if len(args) > 1 and args[1] == "--":
                extra_args = args[2:]
            else:
                extra_args = args[1:]
        else:
            extra_args = args

    # Build pytest command
    cmd = [sys.executable, "-m", "pytest"]

    mode_map = {
        "tier1": ["tests/", "-m", "tier1 and not tier2"],
        "tier2": ["tests/tier2/", "-m", "tier2"],
        "all": ["tests/"],
        "unicode": ["tests/", "-m", "unicode"],
        "unit": ["tests/unit/"],
        "integration": ["tests/integration/"],
    }

    cmd.extend(mode_map.get(mode, mode_map["tier1"]))
    cmd.extend(extra_args)

    print(f"Running: {' '.join(cmd)}")
    return subprocess.call(cmd)


if __name__ == "__main__":
    sys.exit(main())

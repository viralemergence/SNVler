#!/usr/bin/env python3
"""Compatibility wrapper for SNVler primer BED generation."""

import sys

from snvler.cli import main


if __name__ == "__main__":
    raise SystemExit(main(["primers-to-bed", *sys.argv[1:]]))

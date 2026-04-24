#!/usr/bin/env python3
"""Compatibility wrapper for SNVler batch assembly."""

import sys

from snvler.cli import main


if __name__ == "__main__":
    raise SystemExit(main(["batch", *sys.argv[1:]]))

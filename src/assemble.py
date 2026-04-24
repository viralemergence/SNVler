#!/usr/bin/env python3
"""Compatibility wrapper for the packaged SNVler CLI."""

import sys

from snvler.cli import main


if __name__ == "__main__":
    raise SystemExit(main(["assemble", *sys.argv[1:]]))

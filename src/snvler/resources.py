from __future__ import annotations

from contextlib import contextmanager
from importlib.resources import as_file, files
from pathlib import Path
from typing import Iterator


def data_path(*parts: str):
    """Return a traversable package-data path for bundled references or primers."""
    return files("snvler.data").joinpath(*parts)


@contextmanager
def data_file(*parts: str) -> Iterator[Path]:
    """Yield a real filesystem path for bundled package data."""
    with as_file(data_path(*parts)) as path:
        yield path

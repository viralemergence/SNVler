"""Public API for SNVler."""

from .assembly import (
    AssemblyResult,
    Segment,
    SegmentResult,
    assemble_sample,
    build_segments,
    ensure_tools_available,
)
from .batch import BatchResult, assemble_batch
from .primers import generate_bed_from_excel

__all__ = [
    "AssemblyResult",
    "BatchResult",
    "Segment",
    "SegmentResult",
    "assemble_batch",
    "assemble_sample",
    "build_segments",
    "ensure_tools_available",
    "generate_bed_from_excel",
]

__version__ = "0.1.0"

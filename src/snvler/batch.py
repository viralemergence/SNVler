from __future__ import annotations

import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Sequence

import pandas as pd

from .assembly import (
    AssemblyResult,
    PathInput,
    assemble_sample,
    build_segments,
    sample_stem,
)

LogFn = Callable[[str], None]


@dataclass
class BatchResult:
    sample_results: list[AssemblyResult] = field(default_factory=list)
    combined_report: Path | None = None


def read_sample_file(sample_file: PathInput, column: str = "sample") -> list[str]:
    sample_path = Path(sample_file)
    extension = sample_path.suffix.lower()
    if extension == ".csv":
        frame = pd.read_csv(sample_path)
    elif extension in {".xls", ".xlsx"}:
        frame = pd.read_excel(sample_path)
    else:
        raise ValueError("Unsupported sample file type. Use CSV or Excel.")

    if column not in frame.columns:
        raise ValueError(f"Column '{column}' not found in the sample file.")

    return [str(value) for value in frame[column].dropna().tolist()]


def assemble_batch(
    sample_file: PathInput,
    output_dir: PathInput,
    *,
    references: Sequence[PathInput],
    primer_beds: Sequence[PathInput] | None = None,
    column: str = "sample",
    base_path: PathInput | None = None,
    skip_qc: bool = False,
    skip_masking: bool = False,
    stringent: bool = False,
    combined_reference: bool = False,
    threads: int = 1,
    jobs: int = 1,
    use_index: bool = False,
    index_dir: PathInput | None = None,
    report_mode: str = "flagstat",
    dryrun: bool = False,
    log: LogFn = print,
) -> BatchResult:
    sample_names = read_sample_file(sample_file, column)
    segments = build_segments(references, primer_beds)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    base = Path(base_path) if base_path else None

    batch_result = BatchResult()

    def run_one(sample: str) -> AssemblyResult:
        sample_path = base / sample if base else Path(sample)
        sample_output_dir = output_path / sample_stem(sample_path)
        log(f"Running assembly for sample: {sample_path}")
        return assemble_sample(
            [sample_path],
            segments,
            sample_output_dir,
            skip_qc=skip_qc,
            skip_masking=skip_masking,
            stringent=stringent,
            combined_reference=combined_reference,
            threads=threads,
            use_index=use_index,
            index_dir=index_dir,
            report_mode=report_mode,
            dryrun=dryrun,
            log=log,
        )

    if jobs > 1 and len(sample_names) > 1:
        with ThreadPoolExecutor(max_workers=jobs) as executor:
            batch_result.sample_results.extend(executor.map(run_one, sample_names))
    else:
        for sample in sample_names:
            batch_result.sample_results.append(run_one(sample))

    if len(sample_names) > 1:
        if dryrun:
            log("[DRYRUN] Skipping combined mapping report because per-sample reports are not written.")
        else:
            batch_result.combined_report = combine_mapping_reports(sample_names, output_path, log=log)

    return batch_result


def combine_mapping_reports(
    sample_names: Sequence[str],
    output_dir: PathInput,
    *,
    log: LogFn = print,
) -> Path | None:
    output_path = Path(output_dir)
    combined = []

    for sample in sample_names:
        sample_basename = sample_stem(sample)
        out_dir = output_path / sample_basename
        reports_dir = out_dir / "Reports"
        consensus_dir = out_dir / "Consensus"
        mapping_dir = out_dir / "Mapping"

        if not reports_dir.exists():
            log(f"Skipping missing report directory: {reports_dir}")
            continue

        for report_path in sorted(reports_dir.glob("*_mapping_report.csv")):
            frame = pd.read_csv(report_path)
            segment = report_path.name.replace("_mapping_report.csv", "")
            frame["segment"] = segment
            frame["sample"] = sample_basename
            frame["percent_non_N"] = _percent_non_n(consensus_dir / segment / f"{sample_basename}.fa")
            frame["mean_depth"] = _mean_depth(_find_depth_bam(mapping_dir / segment, sample_basename))
            combined.append(frame)

    if not combined:
        return None

    full_frame = pd.concat(combined, ignore_index=True)
    out_csv = output_path / "combined_mapping_report.csv"
    full_frame.to_csv(out_csv, index=False)
    log(f"Combined mapping report written to {out_csv}")
    return out_csv


def _percent_non_n(fasta_file: Path) -> float:
    if not fasta_file.exists():
        return 0.0

    seq_parts: list[str] = []
    with fasta_file.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                seq_parts.append(line.strip().upper())
    sequence = "".join(seq_parts)
    if not sequence:
        return 0.0
    non_n = sum(1 for char in sequence if char != "N")
    return round(non_n / len(sequence) * 100, 2)


def _find_depth_bam(segment_mapping_dir: Path, sample_basename: str) -> Path | None:
    candidates = [
        segment_mapping_dir / f"{sample_basename}_trimmed_mapped.bam",
        segment_mapping_dir / f"{sample_basename}_mapped.bam",
    ]
    candidates.extend(sorted(segment_mapping_dir.glob(f"{sample_basename}*_mapped.bam")))
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _mean_depth(bam_file: Path | None) -> float:
    if bam_file is None:
        return 0.0

    process = subprocess.Popen(
        ["samtools", "depth", "-a", str(bam_file)],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert process.stdout is not None
    depths = [int(line.split()[2]) for line in process.stdout]
    process.wait()
    if not depths:
        return 0.0
    return round(sum(depths) / len(depths), 1)

from __future__ import annotations

import csv
import gzip
import hashlib
import os
import shlex
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Mapping, Sequence

PathInput = str | os.PathLike[str]
LogFn = Callable[[str], None]

EXTERNAL_TOOLS = ("nanoQC", "porechop", "minimap2", "samtools", "ivar")


@dataclass
class Segment:
    """Reference and optional primer BED for one genome segment."""

    name: str
    reference: PathInput
    primer_bed: PathInput | None = None

    def __post_init__(self) -> None:
        self.reference = Path(self.reference)
        if self.primer_bed is not None:
            self.primer_bed = Path(self.primer_bed)

    @classmethod
    def from_reference(
        cls,
        reference: PathInput,
        primer_bed: PathInput | None = None,
        name: str | None = None,
    ) -> "Segment":
        reference_path = Path(reference)
        return cls(
            name=name or sample_stem(reference_path),
            reference=reference_path,
            primer_bed=primer_bed,
        )


@dataclass
class SegmentResult:
    segment: Segment
    mapping_dir: Path
    consensus_dir: Path
    report_csv: Path
    mapped_bams: list[Path] = field(default_factory=list)
    consensus_bams: list[Path] = field(default_factory=list)


@dataclass
class AssemblyResult:
    input_reads: list[Path]
    output_dir: Path
    qc_dir: Path
    processed_reads: list[Path]
    segment_results: list[SegmentResult] = field(default_factory=list)


def sample_stem(path: PathInput) -> str:
    """Return a stable sample-like stem without truncating IDs at underscores."""
    name = Path(path).name
    lower_name = name.lower()
    for suffix in (
        ".fastq.gz",
        ".fq.gz",
        ".fasta.gz",
        ".fa.gz",
        ".fastq",
        ".fq",
        ".fasta",
        ".fa",
        ".bam",
    ):
        if lower_name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def consensus_sample_name(path: PathInput) -> str:
    """Normalize pipeline output names back to the original sample stem."""
    stem = Path(path).with_suffix("").name
    suffixes = (
        "_trimmed_trimmed_sorted",
        "_trimmed_sorted",
        "_trimmed_mapped",
        "_mapped",
        "_trimmed",
    )
    changed = True
    while changed:
        changed = False
        for suffix in suffixes:
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
                changed = True
                break
    return stem


def build_segments(
    references: Sequence[PathInput],
    primer_beds: Sequence[PathInput] | None = None,
) -> list[Segment]:
    """Build segment configs by pairing references and primer BED files by order."""
    if not references:
        raise ValueError("At least one reference FASTA is required.")

    beds = list(primer_beds or [])
    return [
        Segment.from_reference(reference, primer_bed=beds[index] if index < len(beds) else None)
        for index, reference in enumerate(references)
    ]


def first_fasta_record_id(fasta_file: PathInput) -> str:
    """Return the first FASTA record ID, matching samtools/minimap2 contig naming."""
    return first_fasta_record_info(fasta_file)[0]


def first_fasta_record_info(fasta_file: PathInput) -> tuple[str, int]:
    """Return the first FASTA record ID and sequence length."""
    record_id = None
    length = 0
    with Path(fasta_file).open() as handle:
        for line in handle:
            if line.startswith(">"):
                if record_id is not None:
                    break
                record_id = line[1:].strip().split()[0]
                continue
            if record_id is not None:
                length += len(line.strip())
    if record_id is not None:
        return record_id, length
    raise ValueError(f"No FASTA header found in {fasta_file}")


def segment_contig(segment: Segment) -> str:
    return first_fasta_record_id(segment.reference)


def reference_signature(references: Sequence[PathInput]) -> str:
    """Build a short cache signature for reference paths and file metadata."""
    hasher = hashlib.sha1()
    for reference in references:
        path = Path(reference).resolve()
        hasher.update(str(path).encode())
        if path.exists():
            stat = path.stat()
            hasher.update(str(stat.st_size).encode())
            hasher.update(str(int(stat.st_mtime_ns)).encode())
    return hasher.hexdigest()[:16]


def reference_index_path(reference: PathInput, index_dir: PathInput) -> Path:
    reference_path = Path(reference)
    signature = reference_signature([reference_path])
    return Path(index_dir) / f"{reference_path.stem}.{signature}.mmi"


def build_minimap2_index(
    reference: PathInput,
    index_path: PathInput,
    *,
    dryrun: bool = False,
    log: LogFn = print,
) -> Path:
    """Build a minimap2 index if it is absent or stale."""
    reference_path = Path(reference)
    index = Path(index_path)
    index.parent.mkdir(parents=True, exist_ok=True)

    if dryrun and not reference_path.exists():
        run_command(["minimap2", "-d", index, reference_path], dryrun=dryrun, log=log)
        return index

    needs_build = not index.exists() or index.stat().st_mtime < reference_path.stat().st_mtime
    if needs_build:
        run_command(["minimap2", "-d", index, reference_path], dryrun=dryrun, log=log)
    else:
        log(f"Using cached minimap2 index: {index}")
    return index


def mapping_reference(
    reference: PathInput,
    *,
    use_index: bool = False,
    index_dir: PathInput | None = None,
    dryrun: bool = False,
    log: LogFn = print,
) -> Path:
    reference_path = Path(reference)
    if not use_index:
        return reference_path

    cache_dir = Path(index_dir) if index_dir else reference_path.parent / ".snvler_indexes"
    index_path = reference_index_path(reference_path, cache_dir)
    return build_minimap2_index(reference_path, index_path, dryrun=dryrun, log=log)


def combined_reference_path(segments: Sequence[Segment], cache_dir: PathInput) -> Path:
    signature = reference_signature([segment.reference for segment in segments])
    return Path(cache_dir) / f"combined_reference.{signature}.fasta"


def write_combined_reference(
    segments: Sequence[Segment],
    output_path: PathInput,
    *,
    dryrun: bool = False,
    log: LogFn = print,
) -> Path:
    """Write a multi-record reference FASTA preserving original contig IDs."""
    combined_path = Path(output_path)
    combined_path.parent.mkdir(parents=True, exist_ok=True)
    source_paths = [Path(segment.reference) for segment in segments]
    if combined_path.exists() and all(combined_path.stat().st_mtime >= path.stat().st_mtime for path in source_paths):
        log(f"Using cached combined reference: {combined_path}")
        return combined_path

    if dryrun:
        log(f"[DRYRUN] Would write combined reference to {combined_path}")
        return combined_path

    with combined_path.open("w") as out_handle:
        for path in source_paths:
            text = path.read_text()
            out_handle.write(text)
            if not text.endswith("\n"):
                out_handle.write("\n")
    log(f"Combined reference written to {combined_path}")
    return combined_path


def command_to_string(cmd: Sequence[object] | str) -> str:
    if isinstance(cmd, str):
        return cmd
    return shlex.join(str(part) for part in cmd)


def run_command(
    cmd: Sequence[object] | str,
    *,
    dryrun: bool = False,
    log: LogFn = print,
) -> subprocess.CompletedProcess[str] | None:
    """Run one command or print it in dry-run mode."""
    cmd_str = command_to_string(cmd)
    if dryrun:
        log(f"[DRYRUN] {cmd_str}")
        return None

    if isinstance(cmd, str):
        return subprocess.run(cmd, shell=True, check=True, text=True)
    return subprocess.run([str(part) for part in cmd], check=True, text=True)


def ensure_tools_available(
    tools: Sequence[str] = EXTERNAL_TOOLS,
) -> dict[str, str | None]:
    """Return PATH locations for required external tools."""
    return {tool: shutil.which(tool) for tool in tools}


def missing_tools(tools: Sequence[str] = EXTERNAL_TOOLS) -> list[str]:
    return [tool for tool, path in ensure_tools_available(tools).items() if path is None]


def run_fastqc(
    input_reads: Sequence[PathInput],
    output_dir: PathInput,
    *,
    dryrun: bool = False,
    log: LogFn = print,
) -> None:
    """Run nanoQC for each input read file."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    for read_file in input_reads:
        run_command(["nanoQC", "-o", output_path, read_file], dryrun=dryrun, log=log)
        log(f"QC complete for {read_file}")


def run_nanopore_qc(
    input_reads: Sequence[PathInput],
    output_dir: PathInput,
    *,
    dryrun: bool = False,
    log: LogFn = print,
) -> list[Path]:
    """Trim Nanopore adapters with porechop."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    trimmed_files: list[Path] = []

    for read_file in input_reads:
        read_path = Path(read_file)
        out_file = output_path / f"{sample_stem(read_path)}_trimmed.fastq"
        run_command(["porechop", "-i", read_path, "-o", out_file], dryrun=dryrun, log=log)
        log(f"Porechop trimming complete for {read_path}")
        trimmed_files.append(out_file)

    return trimmed_files


def run_mapping(
    input_reads: Sequence[PathInput],
    reference: PathInput,
    mapping_dir: PathInput,
    *,
    threads: int = 1,
    use_index: bool = False,
    index_dir: PathInput | None = None,
    dryrun: bool = False,
    log: LogFn = print,
) -> list[Path]:
    """Map Nanopore reads to a reference and write sorted BAM files."""
    mapping_path = Path(mapping_dir)
    mapping_path.mkdir(parents=True, exist_ok=True)
    reference_path = Path(reference)
    reference_for_mapping = mapping_reference(
        reference_path,
        use_index=use_index,
        index_dir=index_dir,
        dryrun=dryrun,
        log=log,
    )
    mapped_files: list[Path] = []

    for read_file in input_reads:
        read_path = Path(read_file)
        sample_base = sample_stem(read_path)
        bam_out = mapping_path / f"{sample_base}_mapped.bam"
        pipeline = (
            f"minimap2 -t {threads} -ax map-ont {shlex.quote(str(reference_for_mapping))} "
            f"{shlex.quote(str(read_path))} | "
            f"samtools view -@ {threads} -bS - | "
            f"samtools sort -@ {threads} -o {shlex.quote(str(bam_out))} -"
        )
        run_command(["bash", "-o", "pipefail", "-c", pipeline], dryrun=dryrun, log=log)
        log(f"Mapping complete for {read_path} -> {bam_out}")
        mapped_files.append(bam_out)

    return mapped_files


def run_combined_mapping(
    input_reads: Sequence[PathInput],
    segments: Sequence[Segment],
    output_dir: PathInput,
    *,
    threads: int = 1,
    use_index: bool = False,
    index_dir: PathInput | None = None,
    dryrun: bool = False,
    log: LogFn = print,
) -> dict[str, list[Path]]:
    """Map reads once to a combined reference, then split BAMs by segment contig."""
    output_path = Path(output_dir)
    cache_dir = Path(index_dir) if index_dir else output_path / ".snvler_cache"
    combined_ref = write_combined_reference(
        segments,
        combined_reference_path(segments, cache_dir),
        dryrun=dryrun,
        log=log,
    )
    reference_for_mapping = mapping_reference(
        combined_ref,
        use_index=use_index,
        index_dir=cache_dir,
        dryrun=dryrun,
        log=log,
    )
    combined_mapping_dir = output_path / "Mapping" / "_combined"
    combined_mapping_dir.mkdir(parents=True, exist_ok=True)
    mapped_by_segment = {segment.name: [] for segment in segments}

    for read_file in input_reads:
        read_path = Path(read_file)
        sample_base = sample_stem(read_path)
        combined_bam = combined_mapping_dir / f"{sample_base}_combined_mapped.bam"
        pipeline = (
            f"minimap2 -t {threads} -ax map-ont {shlex.quote(str(reference_for_mapping))} "
            f"{shlex.quote(str(read_path))} | "
            f"samtools view -@ {threads} -bS - | "
            f"samtools sort -@ {threads} -o {shlex.quote(str(combined_bam))} -"
        )
        run_command(["bash", "-o", "pipefail", "-c", pipeline], dryrun=dryrun, log=log)
        run_command(["samtools", "index", "-@", threads, combined_bam], dryrun=dryrun, log=log)
        log(f"Combined mapping complete for {read_path} -> {combined_bam}")

        for segment in segments:
            segment_mapping_dir = output_path / "Mapping" / segment.name
            segment_mapping_dir.mkdir(parents=True, exist_ok=True)
            segment_bam = segment_mapping_dir / f"{sample_base}_mapped.bam"
            contig = segment_contig(segment)
            awk_script = (
                'BEGIN{FS=OFS="\\t"} '
                '/^@SQ/ {keep=0; for (i=1; i<=NF; i++) if ($i=="SN:" contig) keep=1; '
                'if (keep) print; next} '
                '/^@/ {print; next} '
                '{print}'
            )
            split_cmd = (
                f"samtools view -@ {threads} -h {shlex.quote(str(combined_bam))} {shlex.quote(contig)} | "
                f"awk -v contig={shlex.quote(contig)} {shlex.quote(awk_script)} | "
                f"samtools view -@ {threads} -b -o {shlex.quote(str(segment_bam))} -"
            )
            run_command(["bash", "-o", "pipefail", "-c", split_cmd], dryrun=dryrun, log=log)
            run_command(["samtools", "index", "-@", threads, segment_bam], dryrun=dryrun, log=log)
            log(f"Split {segment.name} mappings for {read_path} -> {segment_bam}")
            mapped_by_segment[segment.name].append(segment_bam)

    return mapped_by_segment


def primer_masking(
    bam_files: Sequence[PathInput],
    primer_bed: PathInput,
    *,
    threads: int = 1,
    dryrun: bool = False,
    log: LogFn = print,
) -> list[Path]:
    """Trim primer sequence with iVar and return coordinate-sorted BAM files."""
    primer_bed_path = Path(primer_bed)
    trimmed_sorted_bams: list[Path] = []

    for bam in bam_files:
        bam_path = Path(bam)
        sample_base = bam_path.with_suffix("").name
        if sample_base.endswith("_mapped"):
            sample_base = sample_base[: -len("_mapped")]
        out_dir = bam_path.parent
        trimmed_unsorted = out_dir / f"{sample_base}_trimmed.bam"
        output_prefix = trimmed_unsorted.with_suffix("")

        run_command(
            [
                "ivar",
                "trim",
                "-i",
                bam_path,
                "-b",
                primer_bed_path,
                "-p",
                output_prefix,
                "-q",
                "0",
                "-m",
                "10",
                "-s",
                "4",
                "-e",
            ],
            dryrun=dryrun,
            log=log,
        )
        log(f"Primer masking complete for {bam_path} -> {trimmed_unsorted}")

        trimmed_sorted = trimmed_unsorted.with_name(f"{trimmed_unsorted.stem}_sorted.bam")
        run_command(["samtools", "sort", "-@", threads, "-o", trimmed_sorted, trimmed_unsorted], dryrun=dryrun, log=log)
        run_command(["samtools", "index", "-@", threads, trimmed_sorted], dryrun=dryrun, log=log)
        trimmed_sorted_bams.append(trimmed_sorted)

    return trimmed_sorted_bams


def call_consensus(
    bam_files: Sequence[PathInput],
    reference: PathInput,
    consensus_dir: PathInput,
    *,
    stringent: bool = False,
    dryrun: bool = False,
    log: LogFn = print,
) -> list[Path]:
    """Generate iVar consensus FASTA files from sorted BAM files."""
    consensus_path = Path(consensus_dir)
    consensus_path.mkdir(parents=True, exist_ok=True)
    reference_path = Path(reference)
    consensus_fastas: list[Path] = []

    for bam in bam_files:
        bam_path = Path(bam)
        sample_name = consensus_sample_name(bam_path)
        output_prefix = consensus_path / sample_name
        cmd_mpileup = [
            "samtools",
            "mpileup",
            "-d",
            "0",
            "-A",
            "-Q",
            "0",
            "-f",
            reference_path,
            bam_path,
        ]
        if stringent:
            cmd_ivar = ["ivar", "consensus", "-p", output_prefix, "-q", "0", "-t", "0.75", "-m", "10", "-n", "N"]
        else:
            cmd_ivar = ["ivar", "consensus", "-p", output_prefix, "-q", "0", "-t", "0", "-m", "1", "-n", "N"]

        pipeline = f"{command_to_string(cmd_mpileup)} | {command_to_string(cmd_ivar)}"
        run_command(["bash", "-o", "pipefail", "-c", pipeline], dryrun=dryrun, log=log)

        fasta_file = output_prefix.with_suffix(".fa")
        if not dryrun:
            rename_fasta_headers(fasta_file)
        log(f"Consensus generated for {bam_path} with prefix {output_prefix}")
        consensus_fastas.append(fasta_file)

    return consensus_fastas


def rename_fasta_headers(fasta_file: PathInput) -> None:
    """Rename all FASTA headers to the sample name in place."""
    fasta_path = Path(fasta_file)
    sample_name = fasta_path.stem
    lines = fasta_path.read_text().splitlines(keepends=True)
    rewritten = [f">{sample_name}\n" if line.startswith(">") else line for line in lines]
    fasta_path.write_text("".join(rewritten))


def count_fastq_reads(read_files: Sequence[PathInput]) -> dict[str, int]:
    """Count reads in FASTQ/FASTQ.GZ files by sample stem."""
    counts: dict[str, int] = {}
    for read_file in read_files:
        path = Path(read_file)
        opener = gzip.open if path.name.lower().endswith(".gz") else open
        with opener(path, "rt") as handle:
            line_count = sum(1 for _ in handle)
        stem = sample_stem(path)
        count = line_count // 4
        counts[stem] = count
        counts[consensus_sample_name(stem)] = count
    return counts


def fast_mapped_reads(bam_file: PathInput) -> int:
    """Count mapped reads, preferring indexed BAM stats when available."""
    result = subprocess.run(
        ["samtools", "idxstats", str(bam_file)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        result = subprocess.run(
            ["samtools", "view", "-c", "-F", "4", str(bam_file)],
            capture_output=True,
            text=True,
            check=True,
        )
        return int(result.stdout.strip() or 0)

    mapped = 0
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) >= 3 and fields[0] != "*":
            mapped += int(fields[2])
    return mapped


def generate_mapping_report(
    mapped_files: Sequence[PathInput],
    report_csv: PathInput,
    *,
    total_reads_by_sample: Mapping[str, int] | None = None,
    report_mode: str = "flagstat",
    dryrun: bool = False,
    log: LogFn = print,
) -> Path:
    """Generate a CSV report with read mapping counts from samtools flagstat."""
    report_path = Path(report_csv)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    if dryrun:
        log(f"[DRYRUN] Would generate mapping report at {report_path} for files: {list(mapped_files)}")
        return report_path

    with report_path.open("w", newline="") as csvfile:
        fieldnames = ["sample", "mapped_reads", "percentage_mapped"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for bam in mapped_files:
            bam_path = Path(bam)
            sample = consensus_sample_name(bam_path)
            total_reads = total_reads_by_sample.get(sample) if total_reads_by_sample else None

            if report_mode == "idxstats" and total_reads is not None:
                mapped_reads = fast_mapped_reads(bam_path)
            else:
                result = subprocess.run(
                    ["samtools", "flagstat", str(bam_path)],
                    capture_output=True,
                    text=True,
                    check=True,
                )
                mapped_reads = None
                for line in result.stdout.splitlines():
                    if "in total" in line:
                        total_reads = int(line.split()[0])
                    if "mapped (" in line:
                        mapped_reads = int(line.split()[0])
                        break
            if total_reads and mapped_reads is not None:
                percentage = (mapped_reads / total_reads) * 100
            else:
                percentage = 0.0
            writer.writerow(
                {
                    "sample": sample,
                    "mapped_reads": mapped_reads,
                    "percentage_mapped": f"{percentage:.2f}",
                }
            )
            log(f"Report entry for {sample}: {mapped_reads} mapped ({percentage:.2f}%)")

    log(f"Mapping report saved to {report_path}")
    return report_path


def assemble_sample(
    input_reads: Sequence[PathInput],
    segments: Sequence[Segment],
    output_dir: PathInput,
    *,
    skip_qc: bool = False,
    skip_masking: bool = False,
    stringent: bool = False,
    combined_reference: bool = False,
    threads: int = 1,
    use_index: bool = False,
    index_dir: PathInput | None = None,
    report_mode: str = "flagstat",
    dryrun: bool = False,
    log: LogFn = print,
) -> AssemblyResult:
    """Run the assembly workflow for one sample or one set of input read files."""
    if not input_reads:
        raise ValueError("At least one input read file is required.")
    if not segments:
        raise ValueError("At least one segment is required.")

    input_paths = [Path(read_file) for read_file in input_reads]
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    qc_dir = output_path / "QCReports"
    qc_dir.mkdir(parents=True, exist_ok=True)

    if skip_qc:
        log("Skipping quality control...")
        processed_reads = input_paths
    else:
        log("Starting quality control...")
        run_fastqc(input_paths, qc_dir, dryrun=dryrun, log=log)
        processed_reads = run_nanopore_qc(input_paths, qc_dir, dryrun=dryrun, log=log)

    total_reads_by_sample = None
    if report_mode == "idxstats" and not dryrun:
        total_reads_by_sample = count_fastq_reads(processed_reads)

    result = AssemblyResult(
        input_reads=input_paths,
        output_dir=output_path,
        qc_dir=qc_dir,
        processed_reads=processed_reads,
    )

    combined_mapped_files: dict[str, list[Path]] | None = None
    if combined_reference:
        combined_mapped_files = run_combined_mapping(
            processed_reads,
            segments,
            output_path,
            threads=threads,
            use_index=use_index,
            index_dir=index_dir,
            dryrun=dryrun,
            log=log,
        )

    for segment in segments:
        log(f"Processing segment: {segment.name}")
        segment_mapping_dir = output_path / "Mapping" / segment.name
        segment_consensus_dir = output_path / "Consensus" / segment.name
        segment_report = output_path / "Reports" / f"{segment.name}_mapping_report.csv"

        if combined_mapped_files is None:
            mapped_files = run_mapping(
                processed_reads,
                segment.reference,
                segment_mapping_dir,
                threads=threads,
                use_index=use_index,
                index_dir=index_dir,
                dryrun=dryrun,
                log=log,
            )
        else:
            mapped_files = combined_mapped_files[segment.name]

        if skip_masking or segment.primer_bed is None:
            log(f"Skipping masking step for segment {segment.name}")
            consensus_bams = mapped_files
        else:
            consensus_bams = primer_masking(
                mapped_files,
                segment.primer_bed,
                threads=threads,
                dryrun=dryrun,
                log=log,
            )
            log(f"Masking completed for segment {segment.name}")

        call_consensus(
            consensus_bams,
            segment.reference,
            segment_consensus_dir,
            stringent=stringent,
            dryrun=dryrun,
            log=log,
        )
        generate_mapping_report(
            mapped_files,
            segment_report,
            total_reads_by_sample=total_reads_by_sample,
            report_mode=report_mode,
            dryrun=dryrun,
            log=log,
        )

        result.segment_results.append(
            SegmentResult(
                segment=segment,
                mapping_dir=segment_mapping_dir,
                consensus_dir=segment_consensus_dir,
                report_csv=segment_report,
                mapped_bams=mapped_files,
                consensus_bams=consensus_bams,
            )
        )

    return result

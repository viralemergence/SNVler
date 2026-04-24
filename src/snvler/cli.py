from __future__ import annotations

import argparse
import sys

from .assembly import build_segments, ensure_tools_available, assemble_sample
from .batch import assemble_batch
from .primers import generate_bed_from_excel


def main(argv: list[str] | None = None) -> int:
    args_list = list(sys.argv[1:] if argv is None else argv)
    if args_list and args_list[0].startswith("--") and args_list[0] not in {"--help", "-h"}:
        args_list = ["assemble", *args_list]

    parser = build_parser()
    args = parser.parse_args(args_list)
    return args.func(args)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="snvler",
        description="Nanopore viral genome assembly pipeline.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    assemble_parser = subparsers.add_parser("assemble", help="Assemble one sample or read set.")
    add_assembly_args(assemble_parser)
    assemble_parser.set_defaults(func=run_assemble)

    batch_parser = subparsers.add_parser("batch", help="Assemble samples listed in a CSV or Excel file.")
    batch_parser.add_argument("--samplefile", required=True, help="CSV or Excel file containing sample file names.")
    batch_parser.add_argument("--col", default="sample", help="Column containing sample file names.")
    batch_parser.add_argument("--base_path", "--base-path", dest="base_path", help="Directory containing sample files.")
    add_assembly_args(batch_parser, include_input=False)
    batch_parser.set_defaults(func=run_batch)

    primers_parser = subparsers.add_parser("primers-to-bed", help="Generate a primer BED file from Excel.")
    primers_parser.add_argument("--excel_file", "--excel-file", dest="excel_file", required=True)
    primers_parser.add_argument("--fasta_file", "--fasta-file", dest="fasta_file", required=True)
    primers_parser.add_argument("--bed_output", "--bed-output", dest="bed_output", required=True)
    primers_parser.set_defaults(func=run_primers_to_bed)

    doctor_parser = subparsers.add_parser("doctor", help="Check external command-line tool availability.")
    doctor_parser.set_defaults(func=run_doctor)

    return parser


def add_assembly_args(parser: argparse.ArgumentParser, *, include_input: bool = True) -> None:
    if include_input:
        parser.add_argument("--input", nargs="+", required=True, help="Input FASTQ/FAST5 files.")
    parser.add_argument("--references", nargs="+", required=True, help="Reference FASTA files for each segment.")
    parser.add_argument(
        "--primer_bed",
        "--primer-bed",
        dest="primer_bed",
        nargs="+",
        help="Primer BED files matched to references by order.",
    )
    parser.add_argument("--skip_qc", "--skip-qc", dest="skip_qc", action="store_true", help="Skip QC and trimming.")
    parser.add_argument(
        "--skip_masking",
        "--skip-masking",
        dest="skip_masking",
        action="store_true",
        help="Skip primer masking.",
    )
    parser.add_argument("--dryrun", "--dry-run", dest="dryrun", action="store_true", help="Print commands only.")
    parser.add_argument("--output", required=True, help="Output directory.")
    parser.add_argument("--stringent", action="store_true", help="Use stringent consensus parameters.")
    parser.add_argument(
        "--combined-reference",
        action="store_true",
        help="Map once to a combined multi-segment reference, then split BAMs by segment.",
    )
    parser.add_argument("--threads", type=int, default=1, help="Threads passed to minimap2 and samtools.")
    parser.add_argument("--use-index", action="store_true", help="Build and reuse minimap2 .mmi indexes.")
    parser.add_argument("--index-dir", help="Directory for cached minimap2 indexes and combined references.")
    parser.add_argument(
        "--report-mode",
        choices=["flagstat", "idxstats"],
        default="flagstat",
        help="Use flagstat reports or faster indexed/count-based report generation.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Parallel samples for batch mode. Ignored by single-sample assemble.",
    )


def run_assemble(args: argparse.Namespace) -> int:
    segments = build_segments(args.references, args.primer_bed)
    assemble_sample(
        args.input,
        segments,
        args.output,
        skip_qc=args.skip_qc,
        skip_masking=args.skip_masking,
        stringent=args.stringent,
        combined_reference=args.combined_reference,
        threads=args.threads,
        use_index=args.use_index,
        index_dir=args.index_dir,
        report_mode=args.report_mode,
        dryrun=args.dryrun,
    )
    return 0


def run_batch(args: argparse.Namespace) -> int:
    assemble_batch(
        args.samplefile,
        args.output,
        references=args.references,
        primer_beds=args.primer_bed,
        column=args.col,
        base_path=args.base_path,
        skip_qc=args.skip_qc,
        skip_masking=args.skip_masking,
        stringent=args.stringent,
        combined_reference=args.combined_reference,
        threads=args.threads,
        jobs=args.jobs,
        use_index=args.use_index,
        index_dir=args.index_dir,
        report_mode=args.report_mode,
        dryrun=args.dryrun,
    )
    return 0


def run_primers_to_bed(args: argparse.Namespace) -> int:
    generate_bed_from_excel(args.excel_file, args.fasta_file, args.bed_output)
    return 0


def run_doctor(args: argparse.Namespace) -> int:
    status = ensure_tools_available()
    missing = []
    for tool, path in status.items():
        if path is None:
            print(f"{tool}: MISSING")
            missing.append(tool)
        else:
            print(f"{tool}: {path}")

    if missing:
        print("Missing external tools: " + ", ".join(missing))
        return 1
    print("All external tools are available.")
    return 0

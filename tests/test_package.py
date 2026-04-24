from pathlib import Path

import pandas as pd

from snvler import Segment, assemble_sample, generate_bed_from_excel
from snvler.assembly import build_segments, consensus_sample_name, sample_stem
from snvler.cli import main
from snvler.resources import data_file

ROOT = Path(__file__).resolve().parents[1]


def test_sample_name_helpers_preserve_underscores():
    assert sample_stem("sample_A_001.fastq.gz") == "sample_A_001"
    assert consensus_sample_name("sample_A_001_trimmed_mapped.bam") == "sample_A_001"
    assert consensus_sample_name("sample_A_001_trimmed_trimmed_sorted.bam") == "sample_A_001"


def test_build_segments_pairs_references_and_primers():
    segments = build_segments(
        [ROOT / "Reference" / "NM11_SRef.fasta"],
        [ROOT / "Primers" / "SNV_S_primers.bed"],
    )

    assert segments == [
        Segment(
            name="NM11_SRef",
            reference=ROOT / "Reference" / "NM11_SRef.fasta",
            primer_bed=ROOT / "Primers" / "SNV_S_primers.bed",
        )
    ]


def test_assemble_sample_dryrun_uses_repo_reference_and_primer(capsys, tmp_path):
    fastq = ROOT / "tests" / "fixtures" / "sample.fastq"
    segment = Segment.from_reference(
        ROOT / "Reference" / "NM11_SRef.fasta",
        ROOT / "Primers" / "SNV_S_primers.bed",
    )

    result = assemble_sample(
        [fastq],
        [segment],
        tmp_path,
        skip_qc=True,
        dryrun=True,
    )

    captured = capsys.readouterr().out
    assert "minimap2 -t 1 -ax map-ont" in captured
    assert "ivar consensus" in captured
    assert result.segment_results[0].mapped_bams == [
        tmp_path / "Mapping" / "NM11_SRef" / "sample_mapped.bam"
    ]


def test_combined_reference_dryrun_maps_once_for_three_segments(capsys, tmp_path):
    fastq = ROOT / "tests" / "fixtures" / "sample.fastq"
    segments = [
        Segment.from_reference(ROOT / "Reference" / "NM11_LRef.fasta", ROOT / "Primers" / "SNV_L_primers.bed"),
        Segment.from_reference(ROOT / "Reference" / "NM11_MRef.fasta", ROOT / "Primers" / "SNV_M_primers.bed"),
        Segment.from_reference(ROOT / "Reference" / "NM11_SRef.fasta", ROOT / "Primers" / "SNV_S_primers.bed"),
    ]

    assemble_sample(
        [fastq],
        segments,
        tmp_path,
        skip_qc=True,
        skip_masking=True,
        combined_reference=True,
        threads=4,
        use_index=True,
        index_dir=tmp_path / "indexes",
        dryrun=True,
    )

    captured = capsys.readouterr().out
    assert captured.count("minimap2 -t 4 -ax map-ont") == 1
    assert captured.count("samtools view -@ 4 -b -o") == 3
    assert "minimap2 -d" in captured


def test_cli_assemble_dryrun(capsys, tmp_path):
    exit_code = main(
        [
            "assemble",
            "--input",
            str(ROOT / "tests" / "fixtures" / "sample.fastq"),
            "--references",
            str(ROOT / "Reference" / "NM11_SRef.fasta"),
            "--primer-bed",
            str(ROOT / "Primers" / "SNV_S_primers.bed"),
            "--skip-qc",
            "--dry-run",
            "--threads",
            "2",
            "--output",
            str(tmp_path),
        ]
    )

    captured = capsys.readouterr().out
    assert exit_code == 0
    assert "[DRYRUN]" in captured
    assert "minimap2 -t 2" in captured
    assert "sample_mapped.bam" in captured


def test_cli_root_help_shows_subcommands(capsys):
    try:
        main(["--help"])
    except SystemExit as exc:
        assert exc.code == 0

    captured = capsys.readouterr().out
    assert "assemble" in captured
    assert "batch" in captured
    assert "doctor" in captured


def test_generate_bed_from_excel(tmp_path):
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">ref\nAACCGGTT\n")
    primers = tmp_path / "primers.xlsx"
    pd.DataFrame(
        [
            {"Primer_Name": "PAIR_1_L", "Type": "L", "Sequence": "AACC"},
            {"Primer_Name": "PAIR_1_R", "Type": "R", "Sequence": "AACC"},
        ]
    ).to_excel(primers, index=False)

    output = generate_bed_from_excel(primers, fasta, tmp_path / "primers.bed")

    assert output.read_text().splitlines() == [
        "ref\t1\t5\tPAIR_1_LEFT\t60\t+",
        "ref\t5\t9\tPAIR_1_RIGHT\t60\t-",
    ]


def test_bundled_data_paths_are_available():
    with data_file("reference", "NM11_SRef.fasta") as reference:
        assert reference.exists()

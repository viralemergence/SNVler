from __future__ import annotations

from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from .assembly import PathInput


def load_template(fasta_file: PathInput) -> tuple[str, str]:
    """Load the target sequence from a FASTA file."""
    record = next(SeqIO.parse(str(fasta_file), "fasta"))
    return record.id, str(record.seq).upper()


def find_binding_site(template_seq: str, primer_seq: str) -> int:
    """Return the 0-based primer binding start, or -1 if absent."""
    return template_seq.find(primer_seq)


def write_bed_line(
    chrom: str,
    start: int,
    end: int,
    name: str,
    score: int,
    strand: str,
    out_handle,
) -> None:
    out_handle.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")


def generate_bed_from_excel(
    excel_file: PathInput,
    fasta_file: PathInput,
    bed_output: PathInput,
) -> Path:
    """Generate an iVar-compatible primer BED file from an Excel primer table."""
    frame = pd.read_excel(excel_file)
    required_columns = {"Primer_Name", "Type", "Sequence"}
    missing = required_columns.difference(frame.columns)
    if missing:
        raise ValueError(f"Primer Excel file is missing columns: {', '.join(sorted(missing))}")

    chrom, template_seq = load_template(fasta_file)
    output_path = Path(bed_output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w") as bed_out:
        for _, row in frame.iterrows():
            primer_name = str(row["Primer_Name"])
            primer_type = str(row["Type"])
            primer_seq = str(row["Sequence"]).strip().upper()

            if primer_type.upper() == "R":
                query_seq = str(Seq(primer_seq).reverse_complement())
                pos = find_binding_site(template_seq, query_seq)
                strand = "-"
                primer_name = "_".join(primer_name.split("_")[:2]) + "_RIGHT"
            else:
                query_seq = primer_seq
                pos = find_binding_site(template_seq, query_seq)
                strand = "+"
                primer_name = "_".join(primer_name.split("_")[:2]) + "_LEFT"

            if pos == -1:
                print(f"Primer {primer_name} ({primer_type}) not found in the template sequence.")
                continue

            start = pos + 1
            end = pos + len(primer_seq) + 1
            write_bed_line(chrom, start, end, primer_name, 60, strand, bed_out)
            print(f"Primer {primer_name} ({primer_type}) mapped at {start}-{end} on {strand} strand.")

    return output_path

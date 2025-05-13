#!/usr/bin/env python3
#Import packages
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import openpyxl
import argparse
import subprocess
import os
import csv
from pathlib import Path


def load_template(fasta_file):
    """Load the target sequence from a FASTA file."""
    record = next(SeqIO.parse(fasta_file, "fasta"))
    return record.id, str(record.seq)

def find_binding_site(template_seq, primer_seq):
    """Find primer binding site in the template sequence.
       Returns the 0-based start index if found, or -1 if not found."""
    return template_seq.find(primer_seq)

def write_bed_line(chrom, start, end, name, score, strand, out_handle):
    """Write a single BED line to out_handle."""
    out_handle.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")


def generate_bed_from_excel(excel_file, fasta_file, bed_output):
    # Load primer pairs from Excel file.
    # The Excel file is assumed to have columns like:
    # 'Primer_Name', 'Primer_Type' (L or R), and 'Sequence'
    df = pd.read_excel(excel_file)
    
    chrom, template_seq = load_template(fasta_file)
    
    with open(bed_output, "w") as bed_out:
        # Iterate over each row in the Excel file.
        for index, row in df.iterrows():
            primer_name = row['Primer_Name']  # e.g., SNVSSeg_1_L
            primer_type = row['Type']  # 'L' or 'R'
            primer_seq = row['Sequence'].strip().upper()
            
            # For right primers, we search for the reverse complement.
            if primer_type.upper() == 'R':
                primer_seq_rc = str(Seq(primer_seq).reverse_complement())
                pos = find_binding_site(template_seq, primer_seq_rc)
                strand = '-'  # Right primer binds to the negative strand
                # The binding coordinates in BED: start = pos, end = pos + len(primer)
                start = pos + 1 #1-based index for BED
                end = pos + len(primer_seq) + 1 # 1-based index for BED
                #Delete anything after the second underscore
                primer_name = '_'.join(primer_name.split('_')[:2]) + '_RIGHT'  # Rename for BED output
            else:
                pos = find_binding_site(template_seq, primer_seq)
                strand = '+'  # Left primer binds to the positive strand
                start = pos + 1 #1-based index for BED 
                end = pos + len(primer_seq) + 1 # 1-based index for BED
                primer_name = '_'.join(primer_name.split('_')[:2]) + '_LEFT'  # Rename for BED output
            
            if pos == -1:
                print(f"Primer {primer_name} ({primer_type}) not found in the template sequence.")
                continue  # Primer not found in the template sequence.
            else:
                # Write the BED entry. Here we use a fixed score of 60.
                write_bed_line(chrom, start, end, primer_name, 60, strand, bed_out)
                print(f"Primer {primer_name} ({primer_type}) mapped at {start}-{end} on {strand} strand.")

def main():
    parser = argparse.ArgumentParser(description="Generate BED file from primer pairs in Excel.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--excel_file", required=True, help="Path to the Excel file containing primer pairs.")
    parser.add_argument("--fasta_file", required=True, help="Path to the FASTA file of the target sequence.")
    parser.add_argument("--bed_output", required=True, help="Output path for the generated BED file.")
    
    args = parser.parse_args()
    
    generate_bed_from_excel(args.excel_file, args.fasta_file, args.bed_output)

if __name__ == "__main__":
    main()
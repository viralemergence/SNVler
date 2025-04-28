#!/usr/bin/env python3
import os
import sys
import pandas as pd
import argparse
import subprocess
import openpyxl

def main():
    parser = argparse.ArgumentParser(
        description="Wrapper to run assemble.py for each sample listed in a CSV/Excel file."
    )
    parser.add_argument("--samplefile", required=True,
                        help="Path to CSV or Excel file containing sample file names")
    parser.add_argument("--col", default="sample",
                        help="Column name in the file that contains the sample file names (default: 'sample')")
    parser.add_argument("--base_path", required=False,
                        help="Base path where all sample files are located. The file names from the sample file will be joined to this path.")
    parser.add_argument("--output", required=True,
                        help="Base output directory for results")
    parser.add_argument("--references", nargs="+", required=True,
                        help="Reference FASTA files for each segment")
    parser.add_argument("--primer_bed", nargs="+", required=False,
                        help="Primer BED file(s), if needed")
    parser.add_argument("--skip_qc", action="store_true",
                        help="Skip quality control step")
    parser.add_argument("--skip_masking", action="store_true",
                        help="Skip masking step")
    parser.add_argument("--dryrun", action="store_true",
                        help="Dry run: print commands without executing them")
    parser.add_argument("--assemble_script", default="assemble.py",
                        help="Path to the assemble.py script (default: assemble.py)")
    args = parser.parse_args()

    # Read the sample file based on extension.
    ext = os.path.splitext(args.samplefile)[1].lower()
    if ext == ".csv":
        df = pd.read_csv(args.samplefile)
    elif ext in [".xls", ".xlsx"]:
        df = pd.read_excel(args.samplefile)
    else:
        sys.exit("Unsupported file type for sample file. Use CSV or Excel.")

    if args.col not in df.columns:
        sys.exit(f"Column '{args.col}' not found in the sample file.")

    # Get sample file names from the specified column.
    sample_names = df[args.col].dropna().tolist()

    # Iterate through each sample.
    for sample in sample_names:
        # If a base path is provided, join it with the sample file name.
        if args.base_path:
            sample_path = os.path.join(args.base_path, sample)
        else:
            sample_path = sample

        # Create a sample-specific output directory using the sample's basename.
        sample_basename = os.path.basename(sample).split('.')[0]
        sample_output_dir = os.path.join(args.output, sample_basename)
        os.makedirs(sample_output_dir, exist_ok=True)

        # Build the command to call assemble.py.
        cmd = [
            "python3", args.assemble_script,
            "--input", sample_path,
            "--references"
        ] + args.references

        if args.primer_bed:
            cmd += ["--primer_bed"] + args.primer_bed

        if args.skip_qc:
            cmd.append("--skip_qc")
        if args.skip_masking:
            cmd.append("--skip_masking")
        if args.dryrun:
            cmd.append("--dryrun")

        cmd += ["--output", sample_output_dir]

        # Print the command for verification.
        print("Running command for sample:", sample_path)
        print(" ".join(cmd))

        # Execute the command.
        subprocess.run(cmd, check=True)

    if len(sample_names) > 1:
            combined = []
            for sample in sample_names:
                sample_basename = os.path.basename(sample).split(".")[0]
                out_dir         = os.path.join(args.output, sample_basename)
                reports_dir     = os.path.join(out_dir, "Reports")
                consensus_dir   = os.path.join(out_dir, "Consensus")
                mapping_dir     = os.path.join(out_dir, "Mapping")

                # load each segment report
                for fname in os.listdir(reports_dir):
                    if not fname.endswith("_mapping_report.csv"):
                        continue
                    fpath = os.path.join(reports_dir, fname)
                    dfr   = pd.read_csv(fpath)

                    segment = fname.replace("_mapping_report.csv", "")
                    dfr["segment"] = segment
                    dfr["sample"]  = sample_basename

                    # percent non-N
                    fa = os.path.join(consensus_dir, segment, sample_basename + ".fa")
                    seq = []
                    with open(fa) as fh:
                        for L in fh:
                            if not L.startswith(">"):
                                seq.append(L.strip().upper())
                    seq = "".join(seq)
                    nonN = sum(1 for c in seq if c != "N")
                    pct_nonN = (nonN / len(seq) * 100) if seq else 0.0
                    dfr["percent_non_N"] = round(pct_nonN, 2)

                    # mean depth
                    bam = os.path.join(mapping_dir, segment, sample_basename + "_trimmed_mapped.bam")
                    p = subprocess.Popen(
                        ["samtools", "depth", "-a", bam],
                        stdout=subprocess.PIPE, text=True
                    )
                    depths = [int(line.split()[2]) for line in p.stdout]
                    p.wait()
                    mean_depth = (sum(depths) / len(depths)) if depths else 0.0
                    dfr["mean_depth"] = round(mean_depth, 1)

                    combined.append(dfr)

            if combined:
                full_df = pd.concat(combined, ignore_index=True)
                out_csv = os.path.join(args.output, "combined_mapping_report.csv")
                full_df.to_csv(out_csv, index=False)
                print(f"Combined mapping report written to {out_csv}")

if __name__ == "__main__":
    main()
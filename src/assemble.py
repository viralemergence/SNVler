#!/usr/bin/env python3
import argparse
import subprocess
import os
import csv

# Global flag for dry-run mode.
DRYRUN = False

def run_command(cmd):
    """Execute a command (list or string). In dry-run mode, just print it."""
    if isinstance(cmd, list):
        cmd_str = " ".join(cmd)
    else:
        cmd_str = cmd
    if DRYRUN:
        print("[DRYRUN] Command:", cmd_str)
    else:
        subprocess.run(cmd, check=True)

def run_fastqc(input_files, output_dir):
    """Run QC using nanoQC on each input file."""
    for f in input_files:
        cmd = ["nanoQC", "-o", output_dir, f]
        run_command(cmd)
        print(f"QC complete for {f}")

def run_nanopore_qc(input_files, output_dir):
    """run porechop for adapter trimming."""
    trimmed_files = []
    for f in input_files:
        base = os.path.basename(f).replace(".fastq", "")
        out_file = os.path.join(output_dir, f"{base}_trimmed.fastq")
        cmd = ["porechop", "-i", f, "-o", out_file]
        run_command(cmd)
        print(f"Porechop trimming complete for {f}")
        trimmed_files.append(out_file)
    # In dry-run mode, simulate output if needed.
    if DRYRUN and not trimmed_files:
        trimmed_files = [os.path.join(output_dir, os.path.basename(f).replace(".fastq", "") + "_trimmed.fastq") for f in input_files]
    return trimmed_files

def run_mapping(input_files, reference, mapping_dir):
    """
    Map Nanopore reads against the given reference using minimap2.
    The output is piped through samtools to produce a sorted BAM file.
    """
    mapped_files = []
    for f in input_files:
        sample_base = os.path.basename(f).split(".")[0]
        bam_out = os.path.join(mapping_dir, f"{sample_base}_mapped.bam")
        cmd = (
            f"minimap2 -ax map-ont {reference} {f} | "
            f"samtools view -bS - | samtools sort -o {bam_out} -"
        )
        run_command(["bash", "-c", cmd])
        print(f"Mapping complete for {f} -> {bam_out}")
        mapped_files.append(bam_out)
    # In dry-run mode, simulate output if needed.
    if DRYRUN and not mapped_files:
        mapped_files = [os.path.join(mapping_dir, os.path.basename(f).split(".")[0] + "_mapped.bam") for f in input_files]
    return mapped_files

def primer_masking(bam_files, primer_bed):
    """
    Mask reads in the BAM file that overlap with the regions defined in the primer BED file using iVar.
    Returns a list of masked (trimmed) BAM files.
    """
    trimmed_bam = []
    for bam in bam_files:
        sample_base = os.path.basename(bam).replace("_mapped.bam", "")
        trimmed_bam_out = os.path.join(os.path.dirname(bam), f"{sample_base}_trimmed.bam")
        cmd = ["ivar", "trim",
               "-i", bam,
               "-b", primer_bed,
               "-p", trimmed_bam_out,
               "-q", "0",
               "-m", "10",
               "-s", "4"]
        run_command(cmd)
        print(f"Primer masking complete for {bam} -> {trimmed_bam_out}")
        trimmed_bam.append(trimmed_bam_out)
    return trimmed_bam

def call_consensus(bam_files, reference, consensus_dir):
    """
    For each sorted BAM file, generate a consensus sequence.
    
    Calls:
      samtools mpileup -d 0 -A -Q 0 -f <reference> <bam_file> |
      ivar consensus -p <output_prefix> -q 0 -t 0
    """
    os.makedirs(consensus_dir, exist_ok=True)
    
    for bam in bam_files:
        sample_name = os.path.basename(bam).replace("_mapped.bam", "").replace("_sorted.bam", "")
        output_prefix = os.path.join(consensus_dir, sample_name + "_consensus")
        
        cmd_mpileup = [
            "samtools", "mpileup",
            "-d", "0",
            "-A",
            "-Q", "0",
            "-f", reference,
            bam
        ]
        cmd_ivar = [
            "ivar", "consensus",
            "-p", output_prefix,
            "-q", "0",
            "-t", "0"
        ]
        
        if DRYRUN:
            print("[DRYRUN] Would execute:", " ".join(cmd_mpileup), "|", " ".join(cmd_ivar))
        else:
            try:
                p1 = subprocess.Popen(cmd_mpileup, stdout=subprocess.PIPE)
                p2 = subprocess.Popen(cmd_ivar, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                p1.stdout.close() 
                out, err = p2.communicate()
                
                if p2.returncode != 0:
                    print(f"Error generating consensus for {bam}: {err.decode().strip()}")
                else:
                    print(f"Consensus generated for {bam} with prefix {output_prefix}")
            except Exception as e:
                print(f"An exception occurred while processing {bam}: {e}")

def generate_mapping_report(mapped_files, report_csv):
    """Generate a CSV report with the sample name, number and percentage of mapped reads."""
    if DRYRUN:
        print(f"[DRYRUN] Would generate mapping report at {report_csv} for files: {mapped_files}")
        return

    with open(report_csv, "w", newline="") as csvfile:
        fieldnames = ["sample", "mapped_reads", "percentage_mapped"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for bam in mapped_files:
            sample = os.path.basename(bam).split("_")[0]
            result = subprocess.run(["samtools", "flagstat", bam],
                                    capture_output=True, text=True, check=True)
            total_reads = None
            mapped_reads = None
            for line in result.stdout.splitlines():
                if "in total" in line:
                    total_reads = int(line.split()[0])
                if "mapped (" in line:
                    mapped_reads = int(line.split()[0])
                    break
            if total_reads and mapped_reads:
                perc = (mapped_reads / total_reads) * 100
            else:
                perc = 0.0
            writer.writerow({"sample": sample,
                             "mapped_reads": mapped_reads,
                             "percentage_mapped": f"{perc:.2f}"})
            print(f"Report entry for {sample}: {mapped_reads} mapped ({perc:.2f}%)")
    print(f"Mapping report saved to {report_csv}")

def main():
    global DRYRUN

    parser = argparse.ArgumentParser(
        description="Assembler tool for segmented Hantavirus genomes (Nanopore data only)"
    )
    parser.add_argument("--input", nargs="+", required=True,
                        help="Input FASTQ (or FAST5) files")
    parser.add_argument("--references", nargs="+", required=True,
                        help="Reference FASTA files for each segment")
    parser.add_argument("--primer_bed", nargs="+", required=False,
                        help="Primer BED file(s) for each segment (if needed)")
    parser.add_argument("--skip_qc", action="store_true", 
                        help="Skip quality control step")
    parser.add_argument("--skip_masking", action="store_true", 
                        help="Skip masking step")
    parser.add_argument("--dryrun", action="store_true", 
                        help="Dry run: print commands without executing")
    parser.add_argument("--output", required=True,
                        help="Output directory for results")
    args = parser.parse_args()

    DRYRUN = args.dryrun

    os.makedirs(args.output, exist_ok=True)
    qc_dir = os.path.join(args.output, "QCReports")
    os.makedirs(qc_dir, exist_ok=True)

    # Step 2: Quality control (optional)
    if not args.skip_qc:
        print("Starting quality control...")
        run_fastqc(args.input, qc_dir)
        processed_files = run_nanopore_qc(args.input, qc_dir)
    else:
        print("Skipping quality control...")
        processed_files = args.input

    # Process each segment separately using enumeration to match bed files to segments, be sure to provide the refs and the beds in the same order.
    for i, ref in enumerate(args.references):
        segment = os.path.basename(ref).split('.')[0]
        print(f"Processing segment: {segment}")

        # Create subdirectories for this segment.
        segment_mapping_dir = os.path.join(args.output, "Mapping", segment)
        os.makedirs(segment_mapping_dir, exist_ok=True)
        segment_consensus_dir = os.path.join(args.output, "Consensus", segment)
        os.makedirs(segment_consensus_dir, exist_ok=True)
        segment_report = os.path.join(args.output, "Reports", f"{segment}_mapping_report.csv")
        os.makedirs(os.path.dirname(segment_report), exist_ok=True)
        
        # Step 3: Mapping for this segment.
        mapped_files = run_mapping(processed_files, ref, segment_mapping_dir)
        
        # Step 4: Optional masking using primer_masking with the segment-specific BED file.
        if not args.skip_masking and args.primer_bed:
            try:
                bed_file = args.primer_bed[i]
            except IndexError:
                print(f"Warning: No BED file provided for segment {segment}. Skipping masking for this segment.")
                masked_files = mapped_files
            else:
                masked_files = primer_masking(mapped_files, bed_file)
                print("Masking completed for segment", segment)
        else:
            print("Skipping masking step for segment", segment)
            masked_files = mapped_files
        
        # Step 5: Consensus calling for this segment.
        call_consensus(masked_files, ref, segment_consensus_dir)
        
        # Step 6: Generate mapping report for this segment.
        generate_mapping_report(mapped_files, segment_report)
        


if __name__ == "__main__":
    main()
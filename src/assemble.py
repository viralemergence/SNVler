#!/usr/bin/env python3
import argparse
import subprocess
import os
import csv

def run_fastqc(input_files, output_dir):
    """Run QC (here using nanoQC as an example) on each input file."""
    for f in input_files:
        cmd = ["nanoQC", "-o", output_dir, f]
        subprocess.run(cmd, check=True)
        print(f"QC complete for {f}")

def run_nanopore_qc(input_files, output_dir):
    """For Nanopore reads, run porechop for adapter trimming."""
    trimmed_files = []
    for f in input_files:
        base = os.path.basename(f).replace(".fastq", "")
        out_file = os.path.join(output_dir, f"{base}_trimmed.fastq")
        cmd = ["porechop", "-i", f, "-o", out_file]
        subprocess.run(cmd, check=True)
        print(f"Porechop trimming complete for {f}")
        trimmed_files.append(out_file)
    return trimmed_files

def run_mapping(input_files, reference, primer_bed, mapping_dir):
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
        subprocess.run(["bash", "-c", cmd], check=True)
        print(f"Mapping complete for {f} -> {bam_out}")
        mapped_files.append(bam_out)
    return mapped_files

def call_consensus(mapped_files, reference, consensus_dir):
    """
    For each sorted BAM file, generate a consensus sequence.
    
    Calls:
      samtools mpileup -d 0 -A -Q 0 -f <reference> <bam_file> | 
      ivar consensus -p <output_prefix> -q 0 -t 0
    """
    os.makedirs(consensus_dir, exist_ok=True)
    
    for bam in mapped_files:
        # Derive a sample name from the BAM filename.
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
        
        try:
            p1 = subprocess.Popen(cmd_mpileup, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(cmd_ivar, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
            out, err = p2.communicate()
            
            if p2.returncode != 0:
                print(f"Error generating consensus for {bam}: {err.decode().strip()}")
            else:
                print(f"Consensus generated for {bam} with prefix {output_prefix}")
        except Exception as e:
            print(f"An exception occurred while processing {bam}: {e}")

def generate_mapping_report(mapped_files, report_csv):
    """Generate a CSV report with the sample name, number and percentage of mapped reads."""
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

def run_masking(mask_script, reference, maskfile, maskvcf, output_fasta):
    """Call the masking script (using vcf2fastq functionality) to produce a final FASTA."""
    cmd = ["python3", mask_script, reference, maskfile, maskvcf, output_fasta]
    subprocess.run(cmd, check=True)
    print(f"Final assembled FASTA generated at {output_fasta}")

def main():
    parser = argparse.ArgumentParser(
        description="Assembler tool for segmented Hantavirus genomes (Nanopore data only)"
    )
    parser.add_argument("--input", nargs="+", required=True,
                        help="Input FASTQ (or FAST5) files")
    # Now accept multiple references (one per segment)
    parser.add_argument("--references", nargs="+", required=True,
                        help="Reference FASTA files for each segment")
    parser.add_argument("--primer_bed", required=False,
                        help="Primer BED file (if needed)")
    parser.add_argument("--mask_script", required=False,
                        help="Path to the masking script (e.g., mask_vcf.py)")
    parser.add_argument("--maskfile", required=False,
                        help="Mask file for primer scheme")
    parser.add_argument("--maskvcf", required=False,
                        help="VCF file for masking")
    parser.add_argument("--skip_qc", action="store_true", 
                        help="Skip quality control step")
    parser.add_argument("--skip_masking", action="store_true", 
                        help="Skip masking step")
    parser.add_argument("--output", required=True,
                        help="Output directory for results")
    args = parser.parse_args()

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

    # Process each segment separately
    for ref in args.references:
        segment = os.path.basename(ref).split('.')[0]
        print(f"Processing segment: {segment}")

        # Create subdirectories for this segment
        segment_mapping_dir = os.path.join(args.output, "Mapping", segment)
        os.makedirs(segment_mapping_dir, exist_ok=True)
        segment_consensus_dir = os.path.join(args.output, "Consensus", segment)
        os.makedirs(segment_consensus_dir, exist_ok=True)
        segment_report = os.path.join(args.output, "Reports", f"{segment}_mapping_report.csv")
        os.makedirs(os.path.dirname(segment_report), exist_ok=True)
        
        # Step 3: Mapping for this segment
        mapped_files = run_mapping(processed_files, ref, args.primer_bed, segment_mapping_dir)
        
        # Step 4: Consensus calling for this segment
        call_consensus(mapped_files, ref, segment_consensus_dir)
        
        # Step 5: Generate mapping report for this segment
        generate_mapping_report(mapped_files, segment_report)
        
        # Step 6: Optional masking for this segment
        if not args.skip_masking and args.mask_script and args.maskfile and args.maskvcf:
            fasta_out = os.path.join(args.output, "Masked", f"{segment}_assembled.fasta")
            os.makedirs(os.path.dirname(fasta_out), exist_ok=True)
            run_masking(args.mask_script, ref, args.maskfile, args.maskvcf, fasta_out)
        else:
            print("Skipping masking step for segment", segment)

if __name__ == "__main__":
    main()
# SNVler

SNVler is a Python package and command-line workflow for assembling segmented viral genomes from Nanopore reads. It runs QC, adapter trimming, mapping, primer masking, consensus calling, and mapping reports around common bioinformatics tools.

The package can be used from a terminal, a Python script, or a notebook.

## Requirements

SNVler installs the Python dependencies with `pip`:

- Biopython
- openpyxl
- pandas

The assembly workflow also requires these external command-line tools on `PATH`:

- nanoQC
- porechop
- minimap2
- samtools
- ivar

A Conda/Bioconda environment is the recommended way to install the external tools:

```bash
conda create -n snvler python=3.11 -y
conda activate snvler
conda install -c conda-forge -c bioconda nanoqc porechop minimap2 samtools ivar -y
```

## Installation

For local development:

```bash
python -m pip install -e ".[dev]"
```

From GitHub:

```bash
python -m pip install "git+https://github.com/viralemergence/SNVler.git"
```

After installation, check the external tools:

```bash
snvler doctor
```

## Command-Line Usage

Assemble one sample:

```bash
snvler assemble \
  --input sample.fastq \
  --references Reference/NM11_LRef.fasta Reference/NM11_MRef.fasta Reference/NM11_SRef.fasta \
  --primer-bed Primers/SNV_L_primers.bed Primers/SNV_M_primers.bed Primers/SNV_S_primers.bed \
  --threads 4 \
  --output Results/sample
```

Run without executing external tools:

```bash
snvler assemble \
  --input tests/fixtures/sample.fastq \
  --references Reference/NM11_SRef.fasta \
  --primer-bed Primers/SNV_S_primers.bed \
  --skip-qc \
  --dry-run \
  --output Results/dryrun
```

Assemble samples listed in a CSV or Excel file:

```bash
snvler batch \
  --samplefile samples.xlsx \
  --col sample \
  --base-path Reads \
  --references Reference/NM11_LRef.fasta Reference/NM11_MRef.fasta Reference/NM11_SRef.fasta \
  --primer-bed Primers/SNV_L_primers.bed Primers/SNV_M_primers.bed Primers/SNV_S_primers.bed \
  --combined-reference \
  --threads 4 \
  --jobs 2 \
  --use-index \
  --index-dir Results/indexes \
  --report-mode idxstats \
  --output Results
```

For segmented genomes, `--combined-reference` maps each sample once to a multi-segment reference and then splits the BAM by segment. This is usually faster than mapping the same reads separately to each segment. `--threads` is passed to minimap2 and samtools, `--jobs` runs independent samples in parallel during `snvler batch`, and `--use-index` reuses minimap2 `.mmi` indexes.

Generate an iVar primer BED file:

```bash
snvler primers-to-bed \
  --excel-file Primers/Primers_tiling.xlsx \
  --fasta-file Reference/NM11_SRef.fasta \
  --bed-output Primers/SNV_S_primers.generated.bed
```

## Python Usage

```python
from snvler import Segment, assemble_sample

segments = [
    Segment("L", "Reference/NM11_LRef.fasta", "Primers/SNV_L_primers.bed"),
    Segment("M", "Reference/NM11_MRef.fasta", "Primers/SNV_M_primers.bed"),
    Segment("S", "Reference/NM11_SRef.fasta", "Primers/SNV_S_primers.bed"),
]

result = assemble_sample(
    ["sample.fastq"],
    segments,
    "Results/sample",
    stringent=True,
    combined_reference=True,
    threads=4,
    use_index=True,
    index_dir="Results/indexes",
    report_mode="idxstats",
)
```

Bundled SNV reference and primer files are also available through `snvler.resources`:

```python
from snvler.resources import data_file

with data_file("reference", "NM11_SRef.fasta") as reference:
    print(reference)
```

## Development Checks

```bash
python -m compileall -q src tests
python -m pytest
python -m build
```

## Docker

Build from the repository root:

```bash
docker build -f src/dockerfile -t snvler .
```

Run the CLI:

```bash
docker run --rm snvler doctor
```

# Pipeline user guide

## Overview
Snakemake workflow for ivtFOODIE sequencing data processing.

## Quick start

### 1. Configuration
Edit `config.yaml` and set paths correctly:

```yaml
# Use absolute paths (must start with /)
base_dir: "/dshare/xielab/.../codeRefactoring_test"
raw_data_source: "/dshare/xielab/.../your_raw_fastq_directory"
```

### 2. Verify configuration
Run the verification script:

```bash
bash verify_paths.sh
```

If all checks show ✓, the configuration is valid.

### 3. Run the pipeline

#### Option A: Interactive (recommended)
```bash
bash submit_snakemake.sh
```

This script will:
1. Check that Snakemake is installed
2. Run a dry-run
3. Ask for confirmation, then submit jobs

#### Option B: Non-interactive
```bash
bash run_snakemake.sh
```

Submits jobs without confirmation.

## Pipeline flow

```
Raw FASTQ
    ↓
[setup_rawdata]  ← Symlink FASTQs into 00-rawdata_raw
    ↓
00-rawdata_raw/
    ↓
[split_fastq]  Split FASTQ by part
    ↓
00-rawdata/
    ↓
[process_part]  Parallel QC, alignment, extraction
    ↓
01-outputs/
    ↓
[merge_and_analysis]  Merge BAMs, call peaks, deduplicate, conversion ratio
    ↓
02-merge-outputs/
    ↓
[plot_conversion]  Conversion ratio plot
    ↓
converted_ratio.pdf
```

## Output layout

```
{base_dir}/
├── 00-rawdata_raw/           # Symlinks to raw data (created automatically)
│   ├── .setup_done           # Setup completion flag
│   └── {sample}/
│       ├── {sample}_R1.fastq.gz
│       └── {sample}_R2.fastq.gz
├── 00-rawdata/               # Split per-part FASTQs
├── 01-outputs/               # Intermediate outputs
│   ├── fastqc/
│   ├── trim/
│   └── align/
├── 02-merge-outputs/         # Final outputs
│   ├── align/                # Merged BAM
│   ├── bedfiles/             # Peaks and deduplicated BED
│   ├── info_merge/           # Methylation / conversion info
│   └── ratio/                # Conversion ratio files
│       └── {sample}/
│           ├── {sample}_sites.txt
│           ├── {sample}_ratio.igv.tdf
│           └── {sample}_depth.igv.tdf
├── logs/                     # Logs
│   ├── slurm/               # SLURM job logs
│   ├── setup/
│   ├── split/
│   ├── process/
│   ├── merge/
│   └── plot/
└── converted_ratio.pdf       # Final conversion plot
```

## Configuration notes

### Important `config.yaml` fields

```yaml
# Basic settings
species: mouse          # human, mouse, or fly
n_parts: 32            # Number of splits (tune for data size)
base_dir: "/path/to/output"           # Output root (absolute)
raw_data_source: "/path/to/raw/data"  # Raw FASTQ directory (absolute)

# Resources
resources:
  split_fastq:
    threads: 1
    mem_gb: 20
    partition: "compute_plus"
  
  process_part:
    threads: 4
    mem_gb: 20
    partition: "compute_plus"
  
  merge_and_analysis:
    threads: 16
    mem_gb: 64
    partition: "compute_fat"
```

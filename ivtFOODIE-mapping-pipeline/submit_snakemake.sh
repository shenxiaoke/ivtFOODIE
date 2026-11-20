#!/bin/bash

# Snakemake pipeline submission helper
# Run the Bismark workflow on a SLURM cluster

set -e

# Use this script's directory as working directory
WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$WORK_DIR"

# Read base_dir from config.yaml
BASE_DIR=$(grep "^base_dir:" config.yaml | awk '{print $2}' | tr -d '"' | tr -d "'")
# Ensure BASE_DIR is absolute
if [[ ! "$BASE_DIR" = /* ]]; then
    BASE_DIR="$WORK_DIR/$BASE_DIR"
fi

echo "Working directory: $WORK_DIR"
echo "Base directory: $BASE_DIR"
echo ""

# ============= Settings =============
MAX_JOBS=100              # Max concurrent cluster jobs
LATENCY_WAIT=60          # Filesystem latency wait (seconds)
RERUN_INCOMPLETE=true    # Rerun incomplete jobs
KEEP_GOING=true          # Continue other jobs if one fails

# ============= Conda env (optional) =============
# Uncomment to activate conda, e.g. bismark
# source ~/anaconda3/etc/profile.d/conda.sh
# conda activate bismark

# ============= Check Snakemake =============
if ! command -v snakemake &> /dev/null; then
    echo "Error: snakemake not found in PATH"
    echo "Please install snakemake or activate the correct conda environment"
    echo ""
    echo "To install snakemake:"
    echo "  conda activate bismark"
    echo "  conda install -c bioconda -c conda-forge snakemake"
    exit 1
fi

echo "Using snakemake: $(which snakemake)"
echo "Version: $(snakemake --version)"
echo ""

# ============= Log directories =============
mkdir -p "${BASE_DIR}/logs/slurm"
mkdir -p "${BASE_DIR}/logs/split"
mkdir -p "${BASE_DIR}/logs/process"
mkdir -p "${BASE_DIR}/logs/merge"
mkdir -p "${BASE_DIR}/logs/plot"
mkdir -p "${BASE_DIR}/logs/setup"

# ============= Dry-run =============
echo "========================================="
echo "Running dry-run to check workflow..."
echo "========================================="
echo ""

snakemake -n -p --verbose --reason

if [ $? -ne 0 ]; then
    echo ""
    echo "Error: Dry-run failed. Please fix the errors above."
    exit 1
fi

echo ""
echo "========================================="
echo "Dry-run completed successfully!"
echo "========================================="
echo ""

# ============= Confirm submit =============
read -p "Continue with actual run? (y/n) " -n 1 -r
echo ""

if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Pipeline submission cancelled."
    exit 0
fi

# ============= Submit jobs =============
echo ""
echo "========================================="
echo "Submitting pipeline to SLURM cluster..."
echo "========================================="
echo ""

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="${BASE_DIR}/logs/snakemake_${TIMESTAMP}.log"

# Build snakemake command line
SNAKEMAKE_CMD="snakemake \
    --cluster 'sbatch -c {threads} --mem={resources.mem_gb}G -p {resources.partition} -t {resources.time} -o ${BASE_DIR}/logs/slurm/%j_{rule}_{wildcards}.out -e ${BASE_DIR}/logs/slurm/%j_{rule}_{wildcards}.err' \
    --jobs $MAX_JOBS \
    --latency-wait $LATENCY_WAIT \
    --printshellcmds \
    --reason \
    --stats ${BASE_DIR}/logs/snakemake_stats_${TIMESTAMP}.json"

# Optional flags
if [ "$RERUN_INCOMPLETE" = true ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --rerun-incomplete"
fi

if [ "$KEEP_GOING" = true ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --keep-going"
fi

# Run and tee log
echo "Command: $SNAKEMAKE_CMD"
echo ""
echo "Log file: $LOG_FILE"
echo ""

eval $SNAKEMAKE_CMD 2>&1 | tee $LOG_FILE

# ============= Exit status =============
if [ ${PIPESTATUS[0]} -eq 0 ]; then
    echo ""
    echo "========================================="
    echo "Pipeline completed successfully!"
    echo "========================================="
    echo ""
    echo "Check outputs:"
    echo "  - Sites files: ${BASE_DIR}/02-merge-outputs/ratio/*/*_sites.txt"
    echo "  - Peaks: ${BASE_DIR}/02-merge-outputs/bedfiles/*/*_peaks.narrowPeak"
    echo "  - Conversion plot: ${BASE_DIR}/converted_ratio.pdf"
    echo ""
    echo "Generate report:"
    echo "  snakemake --report report.html"
else
    echo ""
    echo "========================================="
    echo "Pipeline failed or incomplete"
    echo "========================================="
    echo ""
    echo "To rerun failed jobs:"
    echo "  bash submit_snakemake.sh"
    echo ""
    echo "To clean failed outputs:"
    echo "  snakemake --delete-all-output <failed_rule>"
    exit 1
fi

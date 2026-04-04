#!/bin/bash
# Run pipeline without interactive confirmation

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

# SLURM log directory under base_dir
mkdir -p "${BASE_DIR}/logs/slurm"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)

snakemake \
    --cluster "sbatch -c {threads} --mem={resources.mem_gb}G -p {resources.partition} -t {resources.time} -o ${BASE_DIR}/logs/slurm/%j_{rule}_{wildcards}.out -e ${BASE_DIR}/logs/slurm/%j_{rule}_{wildcards}.err" \
    --jobs 100 \
    --latency-wait 60 \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds \
    --reason \
    2>&1 | tee "${BASE_DIR}/logs/snakemake_${TIMESTAMP}.log"

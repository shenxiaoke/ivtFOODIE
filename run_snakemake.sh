#!/bin/bash
# 快速运行脚本（无交互确认）

set -e

# 获取脚本所在的绝对路径作为工作目录
WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$WORK_DIR"

# 从config.yaml读取base_dir
BASE_DIR=$(grep "^base_dir:" config.yaml | awk '{print $2}' | tr -d '"' | tr -d "'")
# 确保BASE_DIR是绝对路径
if [[ ! "$BASE_DIR" = /* ]]; then
    BASE_DIR="$WORK_DIR/$BASE_DIR"
fi

# 创建SLURM日志目录（在base_dir下）
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

#!/bin/bash

# Snakemake Pipeline提交脚本
# 用于在SLURM集群上运行bismark pipeline

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

echo "Working directory: $WORK_DIR"
echo "Base directory: $BASE_DIR"
echo ""

# ============= 配置区域 =============
MAX_JOBS=100              # 最大并行作业数
LATENCY_WAIT=60          # 文件系统延迟等待（秒）
RERUN_INCOMPLETE=true    # 重跑未完成的任务
KEEP_GOING=true          # 一个任务失败时继续其他任务

# ============= 激活环境 =============
# 如果使用conda环境，取消下面的注释
# source ~/anaconda3/etc/profile.d/conda.sh
# conda activate bismark

# ============= 检查Snakemake =============
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

# ============= 创建日志目录 =============
mkdir -p "${BASE_DIR}/logs/slurm"
mkdir -p "${BASE_DIR}/logs/split"
mkdir -p "${BASE_DIR}/logs/process"
mkdir -p "${BASE_DIR}/logs/merge"
mkdir -p "${BASE_DIR}/logs/plot"
mkdir -p "${BASE_DIR}/logs/setup"

# ============= Dry-run检查 =============
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

# ============= 确认提交 =============
read -p "Continue with actual run? (y/n) " -n 1 -r
echo ""

if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Pipeline submission cancelled."
    exit 0
fi

# ============= 提交任务 =============
echo ""
echo "========================================="
echo "Submitting pipeline to SLURM cluster..."
echo "========================================="
echo ""

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="${BASE_DIR}/logs/snakemake_${TIMESTAMP}.log"

# 构建snakemake命令
SNAKEMAKE_CMD="snakemake \
    --cluster 'sbatch -c {threads} --mem={resources.mem_gb}G -p {resources.partition} -t {resources.time} -o ${BASE_DIR}/logs/slurm/%j_{rule}_{wildcards}.out -e ${BASE_DIR}/logs/slurm/%j_{rule}_{wildcards}.err' \
    --jobs $MAX_JOBS \
    --latency-wait $LATENCY_WAIT \
    --printshellcmds \
    --reason \
    --stats ${BASE_DIR}/logs/snakemake_stats_${TIMESTAMP}.json"

# 添加可选参数
if [ "$RERUN_INCOMPLETE" = true ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --rerun-incomplete"
fi

if [ "$KEEP_GOING" = true ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --keep-going"
fi

# 执行并记录日志
echo "Command: $SNAKEMAKE_CMD"
echo ""
echo "Log file: $LOG_FILE"
echo ""

eval $SNAKEMAKE_CMD 2>&1 | tee $LOG_FILE

# ============= 检查执行结果 =============
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

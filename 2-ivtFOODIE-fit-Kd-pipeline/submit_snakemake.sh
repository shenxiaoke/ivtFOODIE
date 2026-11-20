#!/bin/bash

# Snakemake Pipeline提交脚本
# 用于在SLURM集群上运行ivtFOODIE postprocessing pipeline
#
# ⚠️ 重要提示：
# - 不要直接运行 "snakemake --cores 8"，这会在本地执行任务，不会提交到集群
# - 必须使用此脚本来提交任务到集群
# - 此脚本会自动使用 --cluster 参数将任务提交到 SLURM 集群
# - 每个config文件使用独立的工作目录，避免锁文件冲突

set -e

# 获取脚本所在的绝对路径作为工作目录
WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$WORK_DIR"

# ============= 参数解析 =============
CONFIG_FILE=""
DRY_RUN=false
PARTITION=""
MAX_JOBS=100
LATENCY_WAIT=60
NON_INTERACTIVE=false
CONDA_ENV_PREFIX="/share/home/shenk/local/anaconda3/envs/bismark"
# 兼容旧调用；若显式指定，只接受 bismark 环境中的 Snakemake。
: "${SNAKEMAKE_BIN:=}"

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        --snakemake)
            SNAKEMAKE_BIN="$2"
            shift 2
            ;;
        -c|--config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        -n|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -p|--partition)
            PARTITION="$2"
            shift 2
            ;;
        -j|--jobs)
            MAX_JOBS="$2"
            shift 2
            ;;
        -y|--yes|--non-interactive)
            NON_INTERACTIVE=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 -c CONFIG_FILE [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  -c, --config FILE     Config file path (required)"
            echo "      --snakemake PATH  Compatibility option; must point into the bismark environment"
            echo "  -n, --dry-run         Run dry-run only (no actual execution)"
            echo "  -p, --partition NAME SLURM partition name (default: from config or 'normal')"
            echo "  -j, --jobs N          Maximum number of parallel jobs (default: 100)"
            echo "  -y, --yes, --non-interactive  Skip interactive confirmation (auto-confirm)"
            echo "  -h, --help            Show this help message"
            echo ""
            echo "Example:"
            echo "  $0 -c config/config-aid.yaml"
            echo "  $0 -c config/config-ssb.yaml -p normal -j 50"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

# 检查必需参数
if [ -z "$CONFIG_FILE" ]; then
    echo "Error: Config file is required"
    echo "Use -h or --help for usage information"
    exit 1
fi

# 检查config文件是否存在
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file not found: $CONFIG_FILE"
    exit 1
fi

# ============= 统一使用 bismark Conda 环境 =============
if [ ! -d "$CONDA_ENV_PREFIX" ]; then
    echo "Error: Required Conda environment was not found:"
    echo "  $CONDA_ENV_PREFIX"
    exit 1
fi

EXPECTED_SNAKEMAKE="$CONDA_ENV_PREFIX/bin/snakemake"

if [ ! -x "$EXPECTED_SNAKEMAKE" ] || [ ! -x "$CONDA_ENV_PREFIX/bin/python" ]; then
    echo "Error: The bismark environment does not contain executable python and snakemake."
    exit 1
fi

# 若调用方显式传入 Snakemake，只允许使用同一个 bismark 环境，避免
# controller 与下游规则落入不同的 Python/工具环境。
if [ -n "$SNAKEMAKE_BIN" ] && [ ! -x "$SNAKEMAKE_BIN" ]; then
    echo "Error: Provided Snakemake is not executable: $SNAKEMAKE_BIN"
    exit 1
fi

if [ -n "$SNAKEMAKE_BIN" ] && [ "$(realpath "$SNAKEMAKE_BIN")" != "$(realpath "$EXPECTED_SNAKEMAKE")" ]; then
    echo "Error: This workflow is configured to use only the bismark Conda environment."
    echo "Expected Snakemake: $EXPECTED_SNAKEMAKE"
    echo "Provided Snakemake: $SNAKEMAKE_BIN"
    exit 1
fi

SNAKEMAKE_BIN="$EXPECTED_SNAKEMAKE"

# sbatch 默认导出当前环境；将 bismark/bin 放在 PATH 最前面后，计算节点
# 中规则里的 python、bgzip、tabix、bedtools 等命令都会解析到该环境。
export PATH="$CONDA_ENV_PREFIX/bin:$PATH"
export CONDA_PREFIX="$CONDA_ENV_PREFIX"
export CONDA_DEFAULT_ENV="bismark"
hash -r

# ============= 从config文件名生成工作目录名 =============
# 例如: config/config-aid.yaml -> .snakemake-aid
CONFIG_BASENAME=$(basename "$CONFIG_FILE" .yaml)
# 移除 "config-" 前缀（如果存在）
if [[ "$CONFIG_BASENAME" == config-* ]]; then
    WORKDIR_NAME="${CONFIG_BASENAME#config-}"
else
    WORKDIR_NAME="$CONFIG_BASENAME"
fi
SNAKEMAKE_WORKDIR=".snakemake-${WORKDIR_NAME}"
CONFIG_FILE="$(realpath "$CONFIG_FILE")"
SNAKEFILE="$WORK_DIR/workflow/Snakefile"

echo "========================================="
echo "Snakemake Cluster Submission"
echo "========================================="
echo "Working directory: $WORK_DIR"
echo "Config file: $CONFIG_FILE"
echo "Snakemake workdir: $SNAKEMAKE_WORKDIR"
echo ""

# ============= 从config读取output_dir =============
OUTPUT_DIR=$(grep "^  output_dir:" "$CONFIG_FILE" 2>/dev/null | awk '{print $2}' | tr -d '"' | tr -d "'" || echo "../results")
# 如果output_dir是相对路径，转换为绝对路径
if [[ ! "$OUTPUT_DIR" = /* ]]; then
    OUTPUT_DIR="$WORK_DIR/$OUTPUT_DIR"
fi

# ============= 从config读取分区和时间 =============
# 尝试从config读取分区
if [ -z "$PARTITION" ]; then
    PARTITION=$(grep -A 2 "^cluster:" "$CONFIG_FILE" 2>/dev/null | grep "partition:" | awk '{print $2}' | tr -d '"' | tr -d "'" || echo "")
fi
# 如果还是没有，使用默认值
if [ -z "$PARTITION" ]; then
    PARTITION="normal"
    echo "Warning: Partition not specified, using default: $PARTITION"
    echo "You can specify partition in config file under 'cluster.partition' or use -p option"
fi

# 读取默认运行时间
DEFAULT_TIME=$(grep "^  runtime:" "$CONFIG_FILE" 2>/dev/null | awk '{print $2}' | tr -d '"' | tr -d "'" || echo "7:00:00")

echo "Output directory: $OUTPUT_DIR"
echo "SLURM partition: $PARTITION"
echo "Default runtime: $DEFAULT_TIME"
echo "Max parallel jobs: $MAX_JOBS"
echo ""

echo "Using snakemake: $SNAKEMAKE_BIN"
echo "Version: $($SNAKEMAKE_BIN --version)"
echo "Unified Conda environment: $CONDA_ENV_PREFIX"
echo "Python: $(command -v python)"
echo ""

# ============= 创建日志目录 =============
mkdir -p "${OUTPUT_DIR}/logs/slurm"
mkdir -p logs

# ============= Dry-run检查 =============
echo "========================================="
echo "Running dry-run to check workflow..."
echo "========================================="
echo ""

# 使用独立的工作目录进行dry-run
"$SNAKEMAKE_BIN" \
    --snakefile "$SNAKEFILE" \
    --configfile "$CONFIG_FILE" \
    --directory "$SNAKEMAKE_WORKDIR" \
    -n -p --verbose --reason

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

# 如果是dry-run模式，退出
if [ "$DRY_RUN" = true ]; then
    echo "Dry-run mode: Exiting without actual submission."
    exit 0
fi

# ============= 确认提交 =============
# 检查是否需要交互式确认
if [ "$NON_INTERACTIVE" = true ]; then
    # 非交互式模式，自动确认
    echo "Running in non-interactive mode, auto-confirming submission..."
elif [ -t 0 ]; then
    # 有终端连接，可以进行交互式确认
    read -p "Continue with actual run? (y/n) " -n 1 -r
    echo ""
    
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Pipeline submission cancelled."
        exit 0
    fi
else
    # 没有终端连接（后台运行），自动确认
    echo "Running in non-interactive mode (no terminal), auto-confirming submission..."
fi

# ============= 提交任务 =============
echo ""
echo "========================================="
echo "Submitting pipeline to SLURM cluster..."
echo "========================================="
echo ""

RUN_TAG="${WORKDIR_NAME}_$(date +%Y%m%d_%H%M%S_%N)_pid$$"
LOG_FILE="${OUTPUT_DIR}/logs/snakemake_${RUN_TAG}.log"

# 构建snakemake命令
# 注意：必须使用 --cluster 参数才能提交到集群，直接使用 --cores 会在本地运行
# 使用 --directory 指定独立的工作目录，避免多个实例的锁文件冲突
SNAKEMAKE_CMD=(
    "$SNAKEMAKE_BIN"
    --snakefile "$SNAKEFILE"
    --configfile "$CONFIG_FILE"
    --directory "$SNAKEMAKE_WORKDIR"
    --cluster "sbatch --export=ALL --cpus-per-task={threads} --mem={resources.mem_mb}M -p {resources.partition} -t {resources.time} -o ${OUTPUT_DIR}/logs/slurm/%j_{rule}_{wildcards}.out -e ${OUTPUT_DIR}/logs/slurm/%j_{rule}_{wildcards}.err"
    --jobs "$MAX_JOBS"
    --default-resources "mem_mb=8000"
    --latency-wait "$LATENCY_WAIT"
    --rerun-incomplete
    --keep-going
    --printshellcmds
    --reason
    --stats "${OUTPUT_DIR}/logs/snakemake_stats_${RUN_TAG}.json"
)

# 执行并记录日志
printf 'Command:'
printf ' %q' "${SNAKEMAKE_CMD[@]}"
printf '\n'
echo ""
echo "Log file: $LOG_FILE"
echo ""
echo "Note: Each job will be submitted to SLURM cluster."
echo "Use 'squeue -u \$USER' to check job status."
echo ""

"${SNAKEMAKE_CMD[@]}" 2>&1 | tee "$LOG_FILE"

# ============= 检查执行结果 =============
if [ ${PIPESTATUS[0]} -eq 0 ]; then
    echo ""
    echo "========================================="
    echo "Pipeline submitted successfully!"
    echo "========================================="
    echo ""
    echo "Monitor jobs with: squeue -u \$USER"
    echo "View logs in: ${OUTPUT_DIR}/logs/slurm/"
    echo ""
else
    echo ""
    echo "========================================="
    echo "Pipeline submission completed with errors"
    echo "========================================="
    echo ""
    echo "Check the log file for details: $LOG_FILE"
    exit 1
fi

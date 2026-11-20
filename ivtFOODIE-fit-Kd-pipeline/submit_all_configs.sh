#!/bin/bash

# 批量将多个 Snakemake workflow controller 提交到 SLURM。
#
# 工作方式：
#   本脚本在本地/login节点运行
#       ↓
#   为每个 config 提交一个独立的 SLURM controller job
#       ↓
#   controller job 在计算节点运行 submit_snakemake.sh
#       ↓
#   submit_snakemake.sh / Snakemake 继续提交和监控下游 SLURM jobs
#
# 特点：
#   1. SSH 断开不会终止 Snakemake controller。
#   2. 所有 config 的 controller job 会快速连续提交，不会等待前一个结束。
#   3. 每个 config 使用独立日志。
#   4. 不再需要使用 nohup、& 或本地后台进程。
#   5. 避免使用 eval，参数传递更加安全。

set -Eeuo pipefail


# ============================================================
# 工作目录
# ============================================================

WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$WORK_DIR"


# ============================================================
# 配置区域
# ============================================================

CONFIG_DIR="config"

# 下游 Snakemake jobs 使用的 partition。
# 该参数会传递给 submit_snakemake.sh。
WORKER_PARTITION="compute_plus"

# 运行 Snakemake controller 的 SLURM partition。
# controller 会一直运行到对应 Snakemake workflow 结束。
CONTROLLER_PARTITION="compute_plus"

# 每个 Snakemake workflow 最多同时提交的下游任务数。
MAX_JOBS=100

# controller job 的资源。
# controller 一般不进行主要计算，但需要构建 DAG、提交和监控任务。
CONTROLLER_CPUS=1
CONTROLLER_MEM="4G"

# 必须足够覆盖整个 Snakemake workflow 的运行时间。
CONTROLLER_TIME="7-00:00:00"

DRY_RUN=false
ASSUME_YES=false

# 可选：显式指定 Snakemake 可执行文件。
SNAKEMAKE_BIN=""


# ============================================================
# 辅助函数
# ============================================================

print_usage() {
    cat <<EOF
Usage:
  $0 [OPTIONS]

Options:
  -d, --config-dir DIR
      Config directory.
      Default: config

  -p, --partition NAME
      Partition used by downstream Snakemake jobs.
      This value is passed to submit_snakemake.sh.
      Default: compute_plus

  --controller-partition NAME
      Partition used by the Snakemake controller jobs.
      Default: compute_plus

  --controller-time TIME
      Time limit for each Snakemake controller job.
      Default: 7-00:00:00

  --controller-mem MEM
      Memory requested by each controller job.
      Default: 4G

  --controller-cpus N
      CPUs requested by each controller job.
      Default: 1

  -j, --jobs N
      Maximum number of parallel downstream jobs per config.
      Default: 100

  --snakemake PATH
      Explicit path to the Snakemake executable.

  -n, --dry-run
      Run submit_snakemake.sh in dry-run mode.
      The dry-run controller itself is still submitted through SLURM.

  -y, --yes
      Skip submission confirmation.

  -b, --background
      Deprecated compatibility option.
      Controller jobs are always submitted through SLURM and therefore
      automatically run independently of the SSH session.

  -h, --help
      Show this help message.

Examples:
  $0

  $0 -y

  $0 -p compute_plus -j 100

  $0 \\
      --controller-partition compute_plus \\
      --controller-time 10-00:00:00 \\
      --controller-mem 8G

  $0 \\
      --snakemake /path/to/conda/env/bin/snakemake

EOF
}


require_value() {
    local option="$1"
    local value="${2:-}"

    if [[ -z "$value" ]]; then
        echo "Error: $option requires a value." >&2
        exit 1
    fi
}


# ============================================================
# 参数解析
# ============================================================

while [[ $# -gt 0 ]]; do
    case "$1" in
        --snakemake)
            require_value "$1" "${2:-}"
            SNAKEMAKE_BIN="$2"
            shift 2
            ;;

        -d|--config-dir)
            require_value "$1" "${2:-}"
            CONFIG_DIR="$2"
            shift 2
            ;;

        -p|--partition)
            require_value "$1" "${2:-}"
            WORKER_PARTITION="$2"
            shift 2
            ;;

        --controller-partition)
            require_value "$1" "${2:-}"
            CONTROLLER_PARTITION="$2"
            shift 2
            ;;

        --controller-time)
            require_value "$1" "${2:-}"
            CONTROLLER_TIME="$2"
            shift 2
            ;;

        --controller-mem)
            require_value "$1" "${2:-}"
            CONTROLLER_MEM="$2"
            shift 2
            ;;

        --controller-cpus)
            require_value "$1" "${2:-}"
            CONTROLLER_CPUS="$2"
            shift 2
            ;;

        -j|--jobs)
            require_value "$1" "${2:-}"
            MAX_JOBS="$2"
            shift 2
            ;;

        -n|--dry-run)
            DRY_RUN=true
            shift
            ;;

        -y|--yes)
            ASSUME_YES=true
            shift
            ;;

        -b|--background)
            # 保留该参数以兼容旧命令。
            # 新版本始终通过 sbatch 提交，因此无需本地后台运行。
            shift
            ;;

        -h|--help)
            print_usage
            exit 0
            ;;

        *)
            echo "Error: Unknown option: $1" >&2
            echo "Use -h or --help for usage information." >&2
            exit 1
            ;;
    esac
done


# ============================================================
# 参数检查
# ============================================================

if ! command -v sbatch >/dev/null 2>&1; then
    echo "Error: sbatch command was not found." >&2
    echo "This script must be run on a machine with access to SLURM." >&2
    exit 1
fi

if [[ ! "$MAX_JOBS" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: --jobs must be a positive integer: $MAX_JOBS" >&2
    exit 1
fi

if [[ ! "$CONTROLLER_CPUS" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: --controller-cpus must be a positive integer." >&2
    exit 1
fi

if [[ ! -d "$CONFIG_DIR" ]]; then
    echo "Error: Config directory not found: $CONFIG_DIR" >&2
    exit 1
fi

CONFIG_DIR="$(cd "$CONFIG_DIR" && pwd)"

SUBMIT_SCRIPT="$WORK_DIR/submit_snakemake.sh"

if [[ ! -f "$SUBMIT_SCRIPT" ]]; then
    echo "Error: submit_snakemake.sh not found:" >&2
    echo "  $SUBMIT_SCRIPT" >&2
    exit 1
fi

if [[ ! -r "$SUBMIT_SCRIPT" ]]; then
    echo "Error: submit_snakemake.sh is not readable:" >&2
    echo "  $SUBMIT_SCRIPT" >&2
    exit 1
fi

if [[ -n "$SNAKEMAKE_BIN" && ! -x "$SNAKEMAKE_BIN" ]]; then
    echo "Error: Specified Snakemake executable is not executable:" >&2
    echo "  $SNAKEMAKE_BIN" >&2
    exit 1
fi


# ============================================================
# 查找 config 文件
# ============================================================

mapfile -t CONFIG_FILES < <(
    find "$CONFIG_DIR" \
        -maxdepth 1 \
        -type f \
        -name "config-*.yaml" \
        -print |
        sort
)

if [[ ${#CONFIG_FILES[@]} -eq 0 ]]; then
    echo "Error: No config files found matching:" >&2
    echo "  $CONFIG_DIR/config-*.yaml" >&2
    exit 1
fi


# ============================================================
# 显示任务信息
# ============================================================

echo "============================================================"
echo "Batch Snakemake Controller Submission"
echo "============================================================"
echo "Working directory:          $WORK_DIR"
echo "Config directory:           $CONFIG_DIR"
echo "Number of configs:          ${#CONFIG_FILES[@]}"
echo "Worker partition:           ${WORKER_PARTITION:-<submit script default>}"
echo "Controller partition:       ${CONTROLLER_PARTITION:-<SLURM default>}"
echo "Maximum jobs per config:    $MAX_JOBS"
echo "Controller CPUs:            $CONTROLLER_CPUS"
echo "Controller memory:          $CONTROLLER_MEM"
echo "Controller time limit:      $CONTROLLER_TIME"
echo "Dry-run:                    $DRY_RUN"

if [[ -n "$SNAKEMAKE_BIN" ]]; then
    echo "Snakemake executable:       $SNAKEMAKE_BIN"
else
    echo "Snakemake executable:       resolved by submit_snakemake.sh"
fi

echo ""
echo "Config files:"

for config_file in "${CONFIG_FILES[@]}"; do
    echo "  - $(basename "$config_file")"
done

echo ""


# ============================================================
# 提交前确认
# ============================================================

if [[ "$ASSUME_YES" == false && "$DRY_RUN" == false ]]; then
    if [[ ! -t 0 ]]; then
        echo "Error: Non-interactive execution requires -y or --yes." >&2
        exit 1
    fi

    read -r -p \
        "Submit ${#CONFIG_FILES[@]} independent controller jobs? (y/n) " \
        REPLY

    if [[ ! "$REPLY" =~ ^[Yy]$ ]]; then
        echo "Batch submission cancelled."
        exit 0
    fi
fi


# ============================================================
# 创建本次运行的日志目录
# ============================================================

RUN_ID="$(date +%Y%m%d_%H%M%S)_pid$$"
LOG_DIR="$WORK_DIR/batch_logs/$RUN_ID"

mkdir -p "$LOG_DIR"

JOB_TABLE="$LOG_DIR/submitted_jobs.tsv"

printf \
    "status\tjob_id\tconfig\tcontroller_script\tstdout_log\tstderr_log\n" \
    > "$JOB_TABLE"

echo ""
echo "Batch log directory:"
echo "  $LOG_DIR"
echo ""


# ============================================================
# 为每个 config 创建并提交 controller job
# ============================================================

SUBMITTED=0
FAILED=0

declare -a JOB_IDS=()

echo "============================================================"
echo "Submitting independent Snakemake controller jobs"
echo "============================================================"
echo ""

for index in "${!CONFIG_FILES[@]}"; do
    config_file="${CONFIG_FILES[$index]}"
    config_name="$(basename "$config_file")"
    config_stem="${config_name%.yaml}"

    # 用于 SLURM job name 和日志名。
    safe_name="$(
        printf '%s' "$config_stem" |
            sed 's/[^A-Za-z0-9_.-]/_/g' |
            cut -c 1-80
    )"

    # index 可以避免不同文件名清理后发生重名。
    task_prefix="$(printf '%03d' "$index")_${safe_name}"

    controller_script="$LOG_DIR/${task_prefix}.controller.sh"
    stdout_log="$LOG_DIR/${task_prefix}.controller.%j.out"
    stderr_log="$LOG_DIR/${task_prefix}.controller.%j.err"

    echo "------------------------------------------------------------"
    echo "Config: $config_name"

    # 为当前 config 生成独立的 controller 脚本。
    #
    # 使用 printf %q 将变量安全地写入脚本，避免路径中空格或特殊字符
    # 造成命令解析错误。
    {
        echo '#!/bin/bash -l'
        echo
        echo 'set -Eeuo pipefail'
        echo

        printf 'WORK_DIR=%q\n' "$WORK_DIR"
        printf 'SUBMIT_SCRIPT=%q\n' "$SUBMIT_SCRIPT"
        printf 'CONFIG_FILE=%q\n' "$config_file"
        printf 'CONFIG_NAME=%q\n' "$config_name"
        printf 'MAX_JOBS=%q\n' "$MAX_JOBS"
        printf 'WORKER_PARTITION=%q\n' "$WORKER_PARTITION"
        printf 'SNAKEMAKE_BIN=%q\n' "$SNAKEMAKE_BIN"
        printf 'DRY_RUN=%q\n' "$DRY_RUN"

        cat <<'CONTROLLER_EOF'

cd "$WORK_DIR"

echo "============================================================"
echo "Snakemake controller started"
echo "============================================================"
echo "Date:              $(date '+%Y-%m-%d %H:%M:%S')"
echo "Hostname:          $(hostname)"
echo "SLURM job ID:      ${SLURM_JOB_ID:-unknown}"
echo "Config:            $CONFIG_FILE"
echo "Worker partition:  ${WORKER_PARTITION:-submit-script-default}"
echo "Maximum jobs:      $MAX_JOBS"
echo "Dry-run:           $DRY_RUN"
echo ""

cmd=(
    bash
    "$SUBMIT_SCRIPT"
    -c "$CONFIG_FILE"
    -j "$MAX_JOBS"
)

if [[ -n "$SNAKEMAKE_BIN" ]]; then
    cmd+=(
        --snakemake "$SNAKEMAKE_BIN"
    )
fi

if [[ -n "$WORKER_PARTITION" ]]; then
    cmd+=(
        -p "$WORKER_PARTITION"
    )
fi

if [[ "$DRY_RUN" == true ]]; then
    cmd+=(
        -n
    )
else
    # controller job 没有交互式终端，因此必须自动确认。
    cmd+=(
        -y
    )
fi

echo "Executing command:"

printf '  '
printf '%q ' "${cmd[@]}"
printf '\n\n'

# 不要在这里使用 & 或 nohup。
#
# exec 会使当前 SLURM controller job 直接成为
# submit_snakemake.sh / Snakemake 进程。
#
# 因此：
#   1. Snakemake 会持续监控下游任务。
#   2. SSH 断开不会影响它。
#   3. controller 的退出状态就是 Snakemake 的退出状态。
exec "${cmd[@]}"

CONTROLLER_EOF
    } > "$controller_script"

    chmod 750 "$controller_script"

    SBATCH_ARGS=(
        --parsable
        --export=ALL
        --job-name="smk_${safe_name}"
        --chdir="$WORK_DIR"
        --cpus-per-task="$CONTROLLER_CPUS"
        --mem="$CONTROLLER_MEM"
        --time="$CONTROLLER_TIME"
        --output="$stdout_log"
        --error="$stderr_log"
    )

    if [[ -n "$CONTROLLER_PARTITION" ]]; then
        SBATCH_ARGS+=(
            --partition="$CONTROLLER_PARTITION"
        )
    fi

    # sbatch 本身只负责把 controller 放入队列。
    # 它不会等待 controller 或 Snakemake workflow 结束。
    if sbatch_output="$(
        sbatch \
            "${SBATCH_ARGS[@]}" \
            "$controller_script" \
            2>&1
    )"; then
        # --parsable 可能返回：
        #   123456
        # 或：
        #   123456;cluster_name
        job_id="${sbatch_output%%;*}"

        JOB_IDS+=("$job_id")

        SUBMITTED=$((SUBMITTED + 1))

        printf \
            "submitted\t%s\t%s\t%s\t%s\t%s\n" \
            "$job_id" \
            "$config_file" \
            "$controller_script" \
            "$stdout_log" \
            "$stderr_log" \
            >> "$JOB_TABLE"

        echo "Submitted controller job: $job_id"
        echo "Controller stdout:        $stdout_log"
        echo "Controller stderr:        $stderr_log"
    else
        FAILED=$((FAILED + 1))

        printf \
            "failed\tNA\t%s\t%s\t%s\t%s\n" \
            "$config_file" \
            "$controller_script" \
            "$stdout_log" \
            "$stderr_log" \
            >> "$JOB_TABLE"

        echo "Failed to submit controller job."
        echo "sbatch error:"
        echo "  $sbatch_output"
    fi

    echo ""
done


# ============================================================
# 总结
# ============================================================

echo "============================================================"
echo "Submission Summary"
echo "============================================================"
echo "Total configs:                  ${#CONFIG_FILES[@]}"
echo "Controller jobs submitted:      $SUBMITTED"
echo "Controller submissions failed:  $FAILED"
echo "Job table:                      $JOB_TABLE"
echo "Log directory:                  $LOG_DIR"
echo ""

if [[ ${#JOB_IDS[@]} -gt 0 ]]; then
    JOB_ID_LIST="$(
        IFS=,
        echo "${JOB_IDS[*]}"
    )"

    echo "Controller job IDs:"
    printf '  %s\n' "${JOB_IDS[@]}"
    echo ""

    echo "Monitor these controller jobs:"
    echo "  squeue -j $JOB_ID_LIST"
    echo ""

    echo "Monitor all controller and downstream jobs:"
    echo "  squeue -u \$USER"
    echo ""

    echo "Inspect controller logs:"
    echo "  tail -f $LOG_DIR/*.controller.*.out"
    echo ""
fi

echo "The local batch submission script has finished."
echo "All successfully submitted controllers now belong to SLURM."
echo "Closing the SSH connection will not terminate them."
echo ""

echo "Each controller starts immediately when SLURM allocates resources."
echo "Different configs do not wait for one another."
echo ""

if [[ "$FAILED" -gt 0 ]]; then
    exit 1
fi


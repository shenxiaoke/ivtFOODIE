#!/bin/bash
# 路径验证脚本
# 用于验证所有配置的路径是否正确

set -e

WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$WORK_DIR"

echo "========================================"
echo "路径配置验证"
echo "========================================"
echo ""
echo "工作目录: $WORK_DIR"
echo ""

# 读取config.yaml中的配置
echo "正在读取 config.yaml ..."
echo ""

BASE_DIR=$(grep "^base_dir:" config.yaml | awk '{print $2}' | tr -d '"' | tr -d "'")
RAW_DATA_SOURCE=$(grep "^raw_data_source:" config.yaml | awk '{print $2}' | tr -d '"' | tr -d "'")
SAMPLE_RENAME_SCRIPT=$(grep "sample_rename_script:" config.yaml | awk '{print $2}' | tr -d '"' | tr -d "'")
SCRIPTS_DIR=$(grep "scripts:" config.yaml | awk '{print $2}' | tr -d '"' | tr -d "'")
BISMARK_ENV=$(grep "bismark_env:" config.yaml | awk '{print $2}' | tr -d '"' | tr -d "'")

echo "配置检查："
echo "----------"

# 检查base_dir
echo -n "1. base_dir: $BASE_DIR ... "
if [[ "$BASE_DIR" = /* ]]; then
    echo "✓ (绝对路径)"
else
    echo "✗ (相对路径，需要改为绝对路径)"
    exit 1
fi

# 检查raw_data_source
echo -n "2. raw_data_source: $RAW_DATA_SOURCE ... "
if [[ "$RAW_DATA_SOURCE" = /* ]]; then
    echo "✓ (绝对路径)"
    if [ -d "$RAW_DATA_SOURCE" ]; then
        echo "   └─ 目录存在 ✓"
        # 检查是否有fastq.gz文件
        FASTQ_COUNT=$(find "$RAW_DATA_SOURCE" -maxdepth 1 -name "*.fastq.gz" -o -name "*.fq.gz" | wc -l)
        echo "   └─ 找到 $FASTQ_COUNT 个 fastq.gz 文件"
    else
        echo "   └─ 警告：目录不存在 ✗"
    fi
else
    echo "✗ (相对路径，需要改为绝对路径)"
    exit 1
fi

# 检查sample_rename_script
echo -n "3. sample_rename_script: $SAMPLE_RENAME_SCRIPT ... "
if [[ "$SAMPLE_RENAME_SCRIPT" = /* ]]; then
    echo "✓ (绝对路径)"
    if [ -f "$SAMPLE_RENAME_SCRIPT" ]; then
        echo "   └─ 文件存在 ✓"
        if [ -x "$SAMPLE_RENAME_SCRIPT" ]; then
            echo "   └─ 可执行 ✓"
        else
            echo "   └─ 警告：不可执行，运行: chmod +x $SAMPLE_RENAME_SCRIPT"
        fi
    else
        echo "   └─ 警告：文件不存在 ✗"
    fi
else
    echo "✗ (相对路径，需要改为绝对路径)"
    exit 1
fi

# 检查scripts目录
echo -n "4. scripts目录: $SCRIPTS_DIR ... "
if [[ "$SCRIPTS_DIR" = /* ]]; then
    echo "✓ (绝对路径)"
    if [ -d "$SCRIPTS_DIR" ]; then
        echo "   └─ 目录存在 ✓"
    else
        echo "   └─ 警告：目录不存在 ✗"
    fi
else
    echo "✗ (相对路径，需要改为绝对路径)"
    exit 1
fi

# 检查bismark环境
echo -n "5. bismark_env: $BISMARK_ENV ... "
if [[ "$BISMARK_ENV" = /* ]]; then
    echo "✓ (绝对路径)"
    if [ -d "$BISMARK_ENV" ]; then
        echo "   └─ 目录存在 ✓"
        # 检查关键工具
        if [ -f "$BISMARK_ENV/bismark" ]; then
            echo "   └─ bismark 可执行文件存在 ✓"
        else
            echo "   └─ 警告：找不到 bismark ✗"
        fi
    else
        echo "   └─ 警告：目录不存在 ✗"
    fi
else
    echo "✗ (相对路径，需要改为绝对路径)"
    exit 1
fi

echo ""
echo "必要脚本检查："
echo "-------------"

# 检查必要的脚本文件
REQUIRED_SCRIPTS=(
    "Snakefile"
    "config.yaml"
    "do_process.sh"
    "mergeBam_callPeak_dedup_callRatio.sh"
    "Sample_rename.new.sh"
    "run_snakemake.sh"
    "submit_snakemake.sh"
)

for script in "${REQUIRED_SCRIPTS[@]}"; do
    echo -n "- $script ... "
    if [ -f "$WORK_DIR/$script" ]; then
        echo "✓"
    else
        echo "✗ (缺失)"
        exit 1
    fi
done

echo ""
echo "样本检测测试："
echo "-------------"

# 检测样本
if [ -d "$RAW_DATA_SOURCE" ]; then
    R1_COUNT=$(find "$RAW_DATA_SOURCE" -maxdepth 1 \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" -o -name "*.R1*.fastq.gz" -o -name "*.R1*.fq.gz" \) 2>/dev/null | wc -l)
    echo "- 找到 $R1_COUNT 个 R1 文件"
    if [ $R1_COUNT -gt 0 ]; then
        echo "  样本列表（前5个）："
        find "$RAW_DATA_SOURCE" -maxdepth 1 \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" -o -name "*.R1*.fastq.gz" -o -name "*.R1*.fq.gz" \) 2>/dev/null | head -5 | while read f; do
            base=$(basename "$f")
            # 提取样本名
            sample=$(echo "$base" | sed 's/_R1.*//' | sed 's/\.R1.*//')
            echo "    └─ $sample"
        done
        if [ $R1_COUNT -gt 5 ]; then
            echo "    └─ ... 还有 $((R1_COUNT - 5)) 个样本"
        fi
    fi
else
    echo "✗ raw_data_source 目录不存在，无法检测样本"
fi

echo ""
echo "工具脚本检查："
echo "-------------"

TOOL_SCRIPTS=(
    "tools/split_string.py"
    "tools/deduplication.py"
    "tools/merge_info.py"
    "tools/read_to_converted_site.py"
    "tools/plot.R"
)

for script in "${TOOL_SCRIPTS[@]}"; do
    echo -n "- $script ... "
    if [ -f "$WORK_DIR/$script" ]; then
        echo "✓"
    else
        echo "✗ (缺失)"
    fi
done

echo ""
echo "========================================"
echo "验证完成！"
echo "========================================"
echo ""
echo "下一步操作："
echo "1. 如果所有路径都正确，运行: bash submit_snakemake.sh"
echo "2. 如果需要修改路径，编辑 config.yaml 文件"
echo "3. 如果原始数据目录不存在，请确认路径或准备数据"
echo ""

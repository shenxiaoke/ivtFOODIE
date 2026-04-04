#!/bin/bash
# 快速测试修复后的 Snakefile

echo "========================================"
echo "快速测试 - Snakemake Pipeline"
echo "========================================"
echo ""

# 检查是否在正确的目录
if [ ! -f "Snakefile" ]; then
    echo "错误：请在 pipeline 目录下运行此脚本"
    exit 1
fi

# 检查 snakemake 是否可用
if ! command -v snakemake &> /dev/null; then
    echo "错误：snakemake 命令不可用"
    echo "请先激活 bismark 环境："
    echo "  conda activate bismark"
    exit 1
fi

echo "1. 检查配置文件..."
if [ ! -f "config.yaml" ]; then
    echo "✗ config.yaml 不存在"
    exit 1
fi
echo "✓ config.yaml 存在"

echo ""
echo "2. 读取配置..."
RAW_DATA_SOURCE=$(grep "^raw_data_source:" config.yaml | awk '{print $2}' | tr -d '"' | tr -d "'")
echo "   raw_data_source: $RAW_DATA_SOURCE"

echo ""
echo "3. 检查原始数据目录..."
if [ -d "$RAW_DATA_SOURCE" ]; then
    echo "✓ 目录存在: $RAW_DATA_SOURCE"
    
    # 统计 R1 文件
    R1_COUNT=$(find "$RAW_DATA_SOURCE" -maxdepth 1 \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" -o -name "*.R1*.fastq.gz" -o -name "*.R1*.fq.gz" \) 2>/dev/null | wc -l)
    echo "✓ 找到 $R1_COUNT 个 R1 文件"
    
    if [ $R1_COUNT -eq 0 ]; then
        echo "✗ 警告：没有找到 R1 文件！"
        echo "  请检查文件命名格式"
        exit 1
    fi
else
    echo "✗ 目录不存在: $RAW_DATA_SOURCE"
    echo "  请检查 config.yaml 中的 raw_data_source 配置"
    exit 1
fi

echo ""
echo "4. 运行 Snakemake dry-run..."
echo "----------------------------------------"

snakemake -n 2>&1 | head -50

echo "----------------------------------------"

# 检查是否成功
if snakemake -n &> /dev/null; then
    echo ""
    echo "========================================"
    echo "✓ 测试通过！"
    echo "========================================"
    echo ""
    echo "下一步："
    echo "  - 查看完整的任务计划：snakemake -n -p"
    echo "  - 提交运行：bash submit_snakemake.sh"
else
    echo ""
    echo "========================================"
    echo "✗ 测试失败"
    echo "========================================"
    echo ""
    echo "请查看上面的错误信息"
    echo "或运行以下命令查看详细错误："
    echo "  snakemake -n -p"
fi

#!/bin/bash

# 必须通过环境变量传入IN和OUT
# 使用方式: IN=/path/to/input OUT=/path/to/output bash Sample_rename.new.sh
if [ -z "$IN" ] || [ -z "$OUT" ]; then
    echo "Error: IN and OUT environment variables must be set"
    echo "Usage: IN=/path/to/input OUT=/path/to/output bash $0"
    exit 1
fi

echo "Input directory: $IN"
echo "Output directory: $OUT"

# 检查输入目录是否存在
if [ ! -d "$IN" ]; then
    echo "Error: Input directory does not exist: $IN"
    exit 1
fi

mkdir -p "$OUT"
pwd=$PWD

for f in "$IN"/*.R1*.gz
do
    base=$(basename "$f")
    prefix=${base%%.R1*}   # 截取掉 _R1 后的部分
    subdir="$OUT/$prefix"

    mkdir -p "$subdir"

    # 建立R1链接
    ln -s "$f" "$subdir/${prefix}_R1.fastq.gz"

    # 找到对应的R2文件并建立链接
    r2file="$IN/${base/.R1/.R2}"
    if [ -f "$r2file" ]; then
        ln -s "$r2file" "$subdir/${prefix}_R2.fastq.gz"
    else
        echo "Warning: no matching R2 for $f"
    fi
done


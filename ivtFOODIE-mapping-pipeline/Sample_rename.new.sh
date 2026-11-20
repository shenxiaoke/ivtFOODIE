#!/bin/bash

# IN and OUT must be passed via environment variables
# Usage: IN=/path/to/input OUT=/path/to/output bash Sample_rename.new.sh
if [ -z "$IN" ] || [ -z "$OUT" ]; then
    echo "Error: IN and OUT environment variables must be set"
    echo "Usage: IN=/path/to/input OUT=/path/to/output bash $0"
    exit 1
fi

echo "Input directory: $IN"
echo "Output directory: $OUT"

# Check input directory exists
if [ ! -d "$IN" ]; then
    echo "Error: Input directory does not exist: $IN"
    exit 1
fi

mkdir -p "$OUT"
pwd=$PWD

for f in "$IN"/*.R1*.gz
do
    base=$(basename "$f")
    prefix=${base%%.R1*}   # Strip from .R1 onward (sample prefix)
    subdir="$OUT/$prefix"

    mkdir -p "$subdir"

    # Symlink R1
    ln -s "$f" "$subdir/${prefix}_R1.fastq.gz"

    # Find matching R2 and symlink
    r2file="$IN/${base/.R1/.R2}"
    if [ -f "$r2file" ]; then
        ln -s "$r2file" "$subdir/${prefix}_R2.fastq.gz"
    else
        echo "Warning: no matching R2 for $f"
    fi
done


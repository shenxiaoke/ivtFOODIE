#!/bin/bash
# Path verification script
# Validates that configured paths are correct

set -e

WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$WORK_DIR"

echo "========================================"
echo "Path configuration verification"
echo "========================================"
echo ""
echo "Working directory: $WORK_DIR"
echo ""

# Read settings from config.yaml
echo "Reading config.yaml ..."
echo ""

BASE_DIR=$(grep "^base_dir:" config.yaml | awk '{print $2}' | tr -d '"' | tr -d "'")
RAW_DATA_SOURCE=$(grep "^raw_data_source:" config.yaml | awk '{print $2}' | tr -d '"' | tr -d "'")
SAMPLE_RENAME_SCRIPT=$(grep "sample_rename_script:" config.yaml | awk '{print $2}' | tr -d '"' | tr -d "'")
SCRIPTS_DIR=$(grep "scripts:" config.yaml | awk '{print $2}' | tr -d '"' | tr -d "'")
BISMARK_ENV=$(grep "bismark_env:" config.yaml | awk '{print $2}' | tr -d '"' | tr -d "'")

echo "Configuration checks:"
echo "----------"

# Check base_dir
echo -n "1. base_dir: $BASE_DIR ... "
if [[ "$BASE_DIR" = /* ]]; then
    echo "✓ (absolute path)"
else
    echo "✗ (relative path; use an absolute path)"
    exit 1
fi

# Check raw_data_source
echo -n "2. raw_data_source: $RAW_DATA_SOURCE ... "
if [[ "$RAW_DATA_SOURCE" = /* ]]; then
    echo "✓ (absolute path)"
    if [ -d "$RAW_DATA_SOURCE" ]; then
        echo "   └─ directory exists ✓"
        # Check for fastq.gz files
        FASTQ_COUNT=$(find "$RAW_DATA_SOURCE" -maxdepth 1 -name "*.fastq.gz" -o -name "*.fq.gz" | wc -l)
        echo "   └─ found $FASTQ_COUNT fastq.gz file(s)"
    else
        echo "   └─ warning: directory does not exist ✗"
    fi
else
    echo "✗ (relative path; use an absolute path)"
    exit 1
fi

# Check sample_rename_script
echo -n "3. sample_rename_script: $SAMPLE_RENAME_SCRIPT ... "
if [[ "$SAMPLE_RENAME_SCRIPT" = /* ]]; then
    echo "✓ (absolute path)"
    if [ -f "$SAMPLE_RENAME_SCRIPT" ]; then
        echo "   └─ file exists ✓"
        if [ -x "$SAMPLE_RENAME_SCRIPT" ]; then
            echo "   └─ executable ✓"
        else
            echo "   └─ warning: not executable; run: chmod +x $SAMPLE_RENAME_SCRIPT"
        fi
    else
        echo "   └─ warning: file does not exist ✗"
    fi
else
    echo "✗ (relative path; use an absolute path)"
    exit 1
fi

# Check scripts directory
echo -n "4. scripts directory: $SCRIPTS_DIR ... "
if [[ "$SCRIPTS_DIR" = /* ]]; then
    echo "✓ (absolute path)"
    if [ -d "$SCRIPTS_DIR" ]; then
        echo "   └─ directory exists ✓"
    else
        echo "   └─ warning: directory does not exist ✗"
    fi
else
    echo "✗ (relative path; use an absolute path)"
    exit 1
fi

# Check bismark environment
echo -n "5. bismark_env: $BISMARK_ENV ... "
if [[ "$BISMARK_ENV" = /* ]]; then
    echo "✓ (absolute path)"
    if [ -d "$BISMARK_ENV" ]; then
        echo "   └─ directory exists ✓"
        # Check key tools
        if [ -f "$BISMARK_ENV/bismark" ]; then
            echo "   └─ bismark executable found ✓"
        else
            echo "   └─ warning: bismark not found ✗"
        fi
    else
        echo "   └─ warning: directory does not exist ✗"
    fi
else
    echo "✗ (relative path; use an absolute path)"
    exit 1
fi

echo ""
echo "Required pipeline scripts:"
echo "-------------"

# Check required pipeline scripts
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
        echo "✗ (missing)"
        exit 1
    fi
done

echo ""
echo "Sample detection test:"
echo "-------------"

# Sample detection
if [ -d "$RAW_DATA_SOURCE" ]; then
    R1_COUNT=$(find "$RAW_DATA_SOURCE" -maxdepth 1 \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" -o -name "*.R1*.fastq.gz" -o -name "*.R1*.fq.gz" \) 2>/dev/null | wc -l)
    echo "- Found $R1_COUNT R1 file(s)"
    if [ $R1_COUNT -gt 0 ]; then
        echo "  Sample names (first 5):"
        find "$RAW_DATA_SOURCE" -maxdepth 1 \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" -o -name "*.R1*.fastq.gz" -o -name "*.R1*.fq.gz" \) 2>/dev/null | head -5 | while read f; do
            base=$(basename "$f")
            # Extract sample name
            sample=$(echo "$base" | sed 's/_R1.*//' | sed 's/\.R1.*//')
            echo "    └─ $sample"
        done
        if [ $R1_COUNT -gt 5 ]; then
            echo "    └─ ... and $((R1_COUNT - 5)) more sample(s)"
        fi
    fi
else
    echo "✗ raw_data_source directory missing; cannot detect samples"
fi

echo ""
echo "Tool scripts:"
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
        echo "✗ (missing)"
    fi
done

echo ""
echo "========================================"
echo "Verification complete."
echo "========================================"
echo ""
echo "Next steps:"
echo "1. If all paths are OK, run: bash submit_snakemake.sh"
echo "2. To change paths, edit config.yaml"
echo "3. If the raw data directory is missing, fix the path or add data"
echo ""

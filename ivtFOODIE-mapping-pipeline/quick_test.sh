#!/bin/bash
# Quick test for the Snakefile workflow

echo "========================================"
echo "Quick test - Snakemake pipeline"
echo "========================================"
echo ""

# Ensure we are in the pipeline directory
if [ ! -f "Snakefile" ]; then
    echo "Error: run this script from the pipeline directory"
    exit 1
fi

# Check snakemake is available
if ! command -v snakemake &> /dev/null; then
    echo "Error: snakemake not found"
    echo "Activate the bismark environment first, e.g.:"
    echo "  conda activate bismark"
    exit 1
fi

echo "1. Checking config file..."
if [ ! -f "config.yaml" ]; then
    echo "✗ config.yaml not found"
    exit 1
fi
echo "✓ config.yaml found"

echo ""
echo "2. Reading configuration..."
RAW_DATA_SOURCE=$(grep "^raw_data_source:" config.yaml | awk '{print $2}' | tr -d '"' | tr -d "'")
echo "   raw_data_source: $RAW_DATA_SOURCE"

echo ""
echo "3. Checking raw data directory..."
if [ -d "$RAW_DATA_SOURCE" ]; then
    echo "✓ Directory exists: $RAW_DATA_SOURCE"
    
    # Count R1 files
    R1_COUNT=$(find "$RAW_DATA_SOURCE" -maxdepth 1 \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" -o -name "*.R1*.fastq.gz" -o -name "*.R1*.fq.gz" \) 2>/dev/null | wc -l)
    echo "✓ Found $R1_COUNT R1 file(s)"
    
    if [ $R1_COUNT -eq 0 ]; then
        echo "✗ Warning: no R1 files found!"
        echo "  Check FASTQ naming (expected *_R1*.fastq.gz, etc.)"
        exit 1
    fi
else
    echo "✗ Directory does not exist: $RAW_DATA_SOURCE"
    echo "  Check raw_data_source in config.yaml"
    exit 1
fi

echo ""
echo "4. Running Snakemake dry-run..."
echo "----------------------------------------"

snakemake -n 2>&1 | head -50

echo "----------------------------------------"

# Check dry-run success
if snakemake -n &> /dev/null; then
    echo ""
    echo "========================================"
    echo "✓ Test passed."
    echo "========================================"
    echo ""
    echo "Next:"
    echo "  - Full DAG: snakemake -n -p"
    echo "  - Submit: bash submit_snakemake.sh"
else
    echo ""
    echo "========================================"
    echo "✗ Test failed"
    echo "========================================"
    echo ""
    echo "See errors above"
    echo "Or run for details:"
    echo "  snakemake -n -p"
fi

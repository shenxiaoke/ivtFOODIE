#!/bin/bash

# created by Shen Ke, 2024/03/22
# modified: read parameters from config.yaml to avoid hard-coding

set -e

# Absolute path to this script's directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="$SCRIPT_DIR/config.yaml"

# Parse config.yaml with Python and export variables
eval $(python3 << EOF
import yaml

with open("$CONFIG_FILE", 'r') as f:
    config = yaml.safe_load(f)

# Path settings
print(f"export BISMARK_ENV='{config['paths']['bismark_env']}'")
print(f"export SCRIPTS_DIR='{config['paths']['scripts']}'")
print(f"export METHYL_EXTRACT='{config['paths']['methyl_extract']}'")

# Reference genomes (exported for bash associative arrays)
for sp, path in config['genomes'].items():
    print(f"export GENOME_{sp.upper()}='{path}'")

# Chromosome sizes
for sp, path in config['chromsizes'].items():
    print(f"export CHROMSIZE_{sp.upper()}='{path}'")

# Bismark options
print(f"export BISMARK_INSERT_SIZE='{config['bismark']['insert_size']}'")
print(f"export BISMARK_QUALITY='{config['bismark']['quality']}'")
print(f"export BISMARK_TRIM_THRESHOLD='{config['bismark']['trim_threshold']}'")
# save_unmapped: keep unmapped reads (default true)
save_unmapped = config['bismark'].get('save_unmapped', True)
print(f"export BISMARK_SAVE_UNMAPPED='{'true' if save_unmapped else 'false'}'")

# Trim Galore options
print(f"export TRIM_QUALITY='{config['trim_galore']['quality']}'")
print(f"export TRIM_STRINGENCY='{config['trim_galore']['stringency']}'")
print(f"export TRIM_ERROR_RATE='{config['trim_galore']['error_rate']}'")

# Cutadapt options
print(f"export CUT_R1_LEFT='{config['cutadapt']['cut_r1_left']}'")
print(f"export CUT_R1_RIGHT='{config['cutadapt']['cut_r1_right']}'")
print(f"export CUT_R2_LEFT='{config['cutadapt']['cut_r2_left']}'")
print(f"export CUT_R2_RIGHT='{config['cutadapt']['cut_r2_right']}'")
EOF
)

# Prepend configured bismark env to PATH
export PATH=$BISMARK_ENV:$PATH

HM=$1
sp=$2
ncpu=$3
sample=$4
sample_splitPart=$5

# Genome paths from config
declare -A genomes
genomes["human"]=$GENOME_HUMAN
genomes["mouse"]=$GENOME_MOUSE
genomes["fly"]=$GENOME_FLY
genomes["lambda"]=$GENOME_LAMBDA

# Chromosome size paths from config
declare -A chromsizes
chromsizes["human"]=$CHROMSIZE_HUMAN
chromsizes["mouse"]=$CHROMSIZE_MOUSE
chromsizes["fly"]=$CHROMSIZE_FLY

# Scripts directory from config
scripts=$SCRIPTS_DIR


function clean {
    IN=$HM/00-rawdata
    OUT=$HM/01-outputs
    
    R1=$IN/$sample_splitPart/$(ls $IN/${sample_splitPart} | grep _R1 | grep q.gz$)
    R2=$IN/$sample_splitPart/$(ls $IN/${sample_splitPart} | grep _R2 | grep q.gz$)
    outQC=$OUT/fastqc/$sample_splitPart
    outTrim=$OUT/trim/$sample_splitPart
    mkdir -p $outQC
    mkdir -p $outTrim
    fastqc $IN/$sample_splitPart/*.gz -t $ncpu -o $outQC
    # trim_galore options from config.yaml
    trim_galore -q $TRIM_QUALITY --stringency $TRIM_STRINGENCY -e $TRIM_ERROR_RATE --paired $R1 $R2 --gzip -o $outTrim -j $ncpu

    # Cutadapt step using cutadapt options from config.yaml
    # Skip cutadapt when all cut lengths are zero
    cutR1=$outTrim/$(ls $outTrim | grep val_1.fq.gz)
    cutR2=$outTrim/$(ls $outTrim | grep val_2.fq.gz)

    if [ "$CUT_R1_LEFT" -eq 0 ] && [ "$CUT_R1_RIGHT" -eq 0 ] && [ "$CUT_R2_LEFT" -eq 0 ] && [ "$CUT_R2_RIGHT" -eq 0 ]; then
        # Skip cutadapt; rename trim_galore outputs only
        echo "[INFO] All Cutadapt cut lengths are 0; skipping cutadapt"
        mv $cutR1 $outTrim/${sample_splitPart}_trim_val_1.fq.gz
        mv $cutR2 $outTrim/${sample_splitPart}_trim_val_2.fq.gz
    else
        # Run cutadapt with config values
        # R1: 5' cut then 3' cut (see --cut order)
        cutadapt --cut=-${CUT_R1_RIGHT} --cut=${CUT_R1_LEFT} $cutR1 | gzip -c > $outTrim/${sample_splitPart}_trim_val_1.fq.gz
        # R2: 5' cut then 3' cut
        cutadapt --cut=${CUT_R2_LEFT} --cut=-${CUT_R2_RIGHT} $cutR2 | gzip -c > $outTrim/${sample_splitPart}_trim_val_2.fq.gz
        rm $cutR1
        rm $cutR2
    fi

}

function align {
    OUT=$HM/01-outputs
    
    genome=${genomes[$sp]}
    outTrim=$OUT/trim/$sample_splitPart
    outAlign=$OUT/align/$sample_splitPart
    mkdir -p $outAlign
    
    R1=$outTrim/$(ls $outTrim | grep val_1.fq.gz)
    R2=$outTrim/$(ls $outTrim | grep val_2.fq.gz)
    
    # Bismark from config (-X insert_size: e.g. 2000 in vivo, 500 in vitro)
    # Add --unmapped when save_unmapped is enabled
    if [ "$BISMARK_SAVE_UNMAPPED" = "true" ]; then
        bismark --genome $genome --non_directional -p $ncpu -1 $R1 -2 $R2 -o $outAlign --temp_dir=$outAlign -X $BISMARK_INSERT_SIZE --unmapped
    else
        bismark --genome $genome --non_directional -p $ncpu -1 $R1 -2 $R2 -o $outAlign --temp_dir=$outAlign -X $BISMARK_INSERT_SIZE
    fi
    
    aln=$outAlign/${sample_splitPart}.bam
    mv $outAlign/*bismark_bt2_pe.bam $aln
    aln_sort=$outAlign/${sample_splitPart}.sorted.bam

    samtools sort -@ $ncpu $aln > $aln_sort
    samtools index $aln_sort
    
    # Unmapped reads outputs
    outUnmapped=$OUT/unmapped/$sample_splitPart
    mkdir -p $outUnmapped
    if [ "$BISMARK_SAVE_UNMAPPED" = "true" ]; then
        # Move unmapped reads into dedicated directory
        mv $outAlign/*_unmapped_reads_1.fq.gz $outUnmapped/${sample_splitPart}_unmapped_R1.fq.gz
        mv $outAlign/*_unmapped_reads_2.fq.gz $outUnmapped/${sample_splitPart}_unmapped_R2.fq.gz
    else
        # Empty placeholder gzip files for Snakemake outputs
        echo "" | gzip > $outUnmapped/${sample_splitPart}_unmapped_R1.fq.gz
        echo "" | gzip > $outUnmapped/${sample_splitPart}_unmapped_R2.fq.gz
    fi
}

function call_conversion {
    IN=$HM/01-outputs
    OUT=$HM/02-merge-outputs

    mkdir -p $OUT/split_align/$sample
    set +e
    ln -s $HM/01-outputs/align/$sample_splitPart/${sample_splitPart}.sorted.bam* $OUT/split_align/$sample/
    set -e

    outSplitAlign=$OUT/split_align/$sample
    outSplitBed=$OUT/split_bedfiles/$sample
    mkdir -p $outSplitBed
    
    #dedupaln=$outSplitAlign/${sample}.deduplicated.bam
    undedupaln=$outSplitAlign/${sample_splitPart}.sorted.bam
    res=$outSplitBed/${sample_splitPart}.reads.undedup.bed
    mkdir -p $OUT/info/$sample
    info=$OUT/info/$sample/${sample_splitPart}.call.o.txt
    # methyl_extract: -t trim_threshold (e.g. 15 in vivo, 0 in vitro), -q quality
    # Script path from config
    python $METHYL_EXTRACT -c 1 -fo $res -a $undedupaln -t $BISMARK_TRIM_THRESHOLD -q $BISMARK_QUALITY --ATAC 'Y' > $info
}


function call_track {
    OUT=$HM/02-merge-outputs
    
    outSplitBed=$OUT/split_bedfiles/$sample
    res=$outSplitBed/${sample_splitPart}.reads.undedup.bed
    track=$outSplitBed/${sample_splitPart}.reads.undedup.summary.bed
    trackScore=$outSplitBed/${sample_splitPart}.reads.undedup.summary_score.bed
    TMPDIR=$outSplitBed
    
    # Use absolute script paths
    python $scripts/summary_bed.py $res abcdefghijklmnop  | sort -k1,1 -k2,2n > $track
    python $scripts/summary_bed_score_trim-test_wo_st_ed.py $res abcdefghijklmnop ABCDEFGHIJKLMNOP | awk '{if($3-$2>2)print$0}' | sort -k1,1 -k2,2n > $trackScore
}


clean

align

call_conversion

call_track

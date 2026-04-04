#!/bin/bash

# created by Shen Ke, 2024/03/22
# modified: 从config.yaml读取参数，避免硬编码

set -e

# 获取脚本所在的绝对路径
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="$SCRIPT_DIR/config.yaml"

# 使用Python解析config.yaml并导出变量
eval $(python3 << EOF
import yaml

with open("$CONFIG_FILE", 'r') as f:
    config = yaml.safe_load(f)

# 路径配置
print(f"export BISMARK_ENV='{config['paths']['bismark_env']}'")
print(f"export SCRIPTS_DIR='{config['paths']['scripts']}'")
print(f"export METHYL_EXTRACT='{config['paths']['methyl_extract']}'")

# 参考基因组 - 导出为关联数组格式
for sp, path in config['genomes'].items():
    print(f"export GENOME_{sp.upper()}='{path}'")

# 染色体大小文件
for sp, path in config['chromsizes'].items():
    print(f"export CHROMSIZE_{sp.upper()}='{path}'")

# Bismark参数
print(f"export BISMARK_INSERT_SIZE='{config['bismark']['insert_size']}'")
print(f"export BISMARK_QUALITY='{config['bismark']['quality']}'")
print(f"export BISMARK_TRIM_THRESHOLD='{config['bismark']['trim_threshold']}'")
# save_unmapped: 是否保留未比对的reads（默认为true）
save_unmapped = config['bismark'].get('save_unmapped', True)
print(f"export BISMARK_SAVE_UNMAPPED='{'true' if save_unmapped else 'false'}'")

# Trim_galore参数
print(f"export TRIM_QUALITY='{config['trim_galore']['quality']}'")
print(f"export TRIM_STRINGENCY='{config['trim_galore']['stringency']}'")
print(f"export TRIM_ERROR_RATE='{config['trim_galore']['error_rate']}'")

# Cutadapt参数
print(f"export CUT_R1_LEFT='{config['cutadapt']['cut_r1_left']}'")
print(f"export CUT_R1_RIGHT='{config['cutadapt']['cut_r1_right']}'")
print(f"export CUT_R2_LEFT='{config['cutadapt']['cut_r2_left']}'")
print(f"export CUT_R2_RIGHT='{config['cutadapt']['cut_r2_right']}'")
EOF
)

# 设置PATH使用配置中的bismark环境
export PATH=$BISMARK_ENV:$PATH

HM=$1
sp=$2
ncpu=$3
sample=$4
sample_splitPart=$5

# 使用从config读取的基因组路径
declare -A genomes
genomes["human"]=$GENOME_HUMAN
genomes["mouse"]=$GENOME_MOUSE
genomes["fly"]=$GENOME_FLY
genomes["lambda"]=$GENOME_LAMBDA

# 使用从config读取的染色体大小文件路径
declare -A chromsizes
chromsizes["human"]=$CHROMSIZE_HUMAN
chromsizes["mouse"]=$CHROMSIZE_MOUSE
chromsizes["fly"]=$CHROMSIZE_FLY

# 使用从config读取的scripts路径
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
    # 使用config.yaml中的trim_galore参数
    trim_galore -q $TRIM_QUALITY --stringency $TRIM_STRINGENCY -e $TRIM_ERROR_RATE --paired $R1 $R2 --gzip -o $outTrim -j $ncpu

    # Cutadapt处理 - 使用config.yaml中的cutadapt参数
    # 如果所有cut参数都是0，则跳过cutadapt步骤
    cutR1=$outTrim/$(ls $outTrim | grep val_1.fq.gz)
    cutR2=$outTrim/$(ls $outTrim | grep val_2.fq.gz)

    if [ "$CUT_R1_LEFT" -eq 0 ] && [ "$CUT_R1_RIGHT" -eq 0 ] && [ "$CUT_R2_LEFT" -eq 0 ] && [ "$CUT_R2_RIGHT" -eq 0 ]; then
        # 跳过cutadapt，直接重命名文件
        echo "[INFO] Cutadapt参数全为0，跳过cutadapt步骤"
        mv $cutR1 $outTrim/${sample_splitPart}_trim_val_1.fq.gz
        mv $cutR2 $outTrim/${sample_splitPart}_trim_val_2.fq.gz
    else
        # 执行cutadapt，使用config中的参数
        # R1: --cut=左端 --cut=-右端
        cutadapt --cut=-${CUT_R1_RIGHT} --cut=${CUT_R1_LEFT} $cutR1 | gzip -c > $outTrim/${sample_splitPart}_trim_val_1.fq.gz
        # R2: --cut=左端 --cut=-右端
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
    
    # 使用config.yaml中的bismark参数（-X insert_size: in vivo数据用2000，in vitro用500）
    # 根据配置决定是否添加 --unmapped 参数保留未比对的reads
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
    
    # 处理unmapped reads输出
    outUnmapped=$OUT/unmapped/$sample_splitPart
    mkdir -p $outUnmapped
    if [ "$BISMARK_SAVE_UNMAPPED" = "true" ]; then
        # 将unmapped reads移动到专门的目录
        mv $outAlign/*_unmapped_reads_1.fq.gz $outUnmapped/${sample_splitPart}_unmapped_R1.fq.gz
        mv $outAlign/*_unmapped_reads_2.fq.gz $outUnmapped/${sample_splitPart}_unmapped_R2.fq.gz
    else
        # 创建空的占位文件以满足Snakemake输出要求
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
    # 使用config.yaml中的参数（-t trim_threshold: in vivo用15, in vitro用0; -q quality）
    # 使用config中配置的methyl_extract脚本路径
    python $METHYL_EXTRACT -c 1 -fo $res -a $undedupaln -t $BISMARK_TRIM_THRESHOLD -q $BISMARK_QUALITY --ATAC 'Y' > $info
}


function call_track {
    OUT=$HM/02-merge-outputs
    
    outSplitBed=$OUT/split_bedfiles/$sample
    res=$outSplitBed/${sample_splitPart}.reads.undedup.bed
    track=$outSplitBed/${sample_splitPart}.reads.undedup.summary.bed
    trackScore=$outSplitBed/${sample_splitPart}.reads.undedup.summary_score.bed
    TMPDIR=$outSplitBed
    
    # 使用绝对路径
    python $scripts/summary_bed.py $res abcdefghijklmnop  | sort -k1,1 -k2,2n > $track
    python $scripts/summary_bed_score_trim-test_wo_st_ed.py $res abcdefghijklmnop ABCDEFGHIJKLMNOP | awk '{if($3-$2>2)print$0}' | sort -k1,1 -k2,2n > $trackScore
}


clean

align

call_conversion

call_track

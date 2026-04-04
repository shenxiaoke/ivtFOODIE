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
print(f"export SORT_TMP='{config['paths']['sort_tmp']}'")
print(f"export UCSC_TOOLS='{config['paths']['ucsc_tools']}'")

# 染色体大小文件
for sp, path in config['chromsizes'].items():
    print(f"export CHROMSIZE_{sp.upper()}='{path}'")

# MACS2 gsize参数
for sp, gs in config['gsize'].items():
    print(f"export GSIZE_{sp.upper()}='{gs}'")

# IGV genome参数
for sp, igv in config['igv_genome'].items():
    print(f"export IGV_{sp.upper()}='{igv}'")
EOF
)

export PATH=$BISMARK_ENV:$PATH

HM=$1
sp=$2
ncpu=$3
sample=$4

# 使用从config读取的染色体大小文件路径
declare -A chromsizes
chromsizes["human"]=$CHROMSIZE_HUMAN
chromsizes["mouse"]=$CHROMSIZE_MOUSE
chromsizes["fly"]=$CHROMSIZE_FLY

# 使用从config读取的MACS2 gsize参数
declare -A gsize_array
gsize_array["human"]=$GSIZE_HUMAN
gsize_array["mouse"]=$GSIZE_MOUSE
gsize_array["fly"]=$GSIZE_FLY

# 使用从config读取的IGV genome参数
declare -A igv_tdf_array
igv_tdf_array["human"]=$IGV_HUMAN
igv_tdf_array["mouse"]=$IGV_MOUSE
igv_tdf_array["fly"]=$IGV_FLY

function merge_bam_without_dedup {
    IN=$HM/02-merge-outputs/split_align/$sample
    OUT=$HM/02-merge-outputs
    
    mkdir -p $OUT/align/$sample
    aln=$OUT/align/$sample/${sample}.bam
    aln_sortByName=$OUT/align/$sample/${sample}.sort_name.bam
    
    samtools merge -@ $ncpu $aln $(find $IN -name ${sample}-part*.bam)
    samtools sort -n -T $SORT_TMP -@ $ncpu $aln -o $aln_sortByName
}

function call_peak {
    OUT=$HM/02-merge-outputs

    outAlign=$OUT/align/$sample
    aln=$outAlign/${sample}.sort_name.bam
    outBed=$OUT/bedfiles/$sample
    # frag=$outBed/${sample}_frag.bed
    # track=$outBed/${sample}_coverage.bigWig
    chrom=${chromsizes[$sp]}

    mkdir -p $outBed
    
    peakdir=$outBed
    sample=$sample
    aln=$aln
    chromSize=$chrom
    gsize=${gsize_array[$sp]} # hs for human while mm for mouse 

    mkdir -p ${peakdir}/${sample}
    bedtools bamtobed -i $aln -bedpe | awk 'BEGIN {OFS="\t"} length($1) <= 5 && $8 >= 30 && ($6-$2) >=10 {print $1, $2 + 4, $6 - 5}' | sort -k 1,1 -k 2,2n > ${peakdir}/${sample}/frag.bed

    bedtools genomecov -i ${peakdir}/${sample}/frag.bed -g $chromSize -bg > ${peakdir}/${sample}/coverage.bg
    $UCSC_TOOLS/bedGraphToBigWig ${peakdir}/${sample}/coverage.bg $chromSize ${peakdir}/${sample}/coverage.bw

    macs2 callpeak --nomodel --nolambda --call-summits --keep-dup all -g $gsize -t ${peakdir}/${sample}/frag.bed -f BEDPE -n $sample --outdir ${peakdir}/${sample}
    awk 'BEGIN {OFS="\t"} {if ($5 > 50) print $1, $2, $3, $1"_"$2"_"$3}' ${peakdir}/${sample}/${sample}_peaks.narrowPeak | sort -k 1,1 -k 2,2n | uniq > ${peakdir}/${sample}/high_score_peaks.bed
}

function merge_bed_and_dedup {
    IN=$HM/02-merge-outputs/split_bedfiles
    OUT=$HM/02-merge-outputs/bedfiles
    outBed=$OUT/$sample
    mkdir -p $outBed

    mkdir -p $OUT
    cat $IN/$sample/${sample}-part*.reads.undedup.summary_score.bed > $OUT/$sample/$sample.reads.undedup.summary_score.unsort.bed
    sort -T $SORT_TMP -k1,1 -k2,2n -k3,3n $OUT/$sample/$sample.reads.undedup.summary_score.unsort.bed > $OUT/$sample/$sample.reads.undedup.summary_score.bed

    bed_inf=$OUT/$sample/${sample}.reads.undedup.summary_score.bed
    bed_outf=$OUT/$sample/${sample}.reads.undedup.summary_score.deduplicated.bed
    echo $bed_inf $bed_outf|python $SCRIPT_DIR/tools/deduplication.py
}

function merge_bed_and_dedup_readsStartEnd {
    IN=$HM/02-merge-outputs/split_bedfiles
    OUT=$HM/02-merge-outputs/bedfiles
    outBed=$OUT/$sample
    mkdir -p $outBed

    mkdir -p $OUT
    cat $IN/$sample/${sample}-part*.reads.undedup.summary_score.bed > $OUT/$sample/$sample.reads.undedup.summary_score.unsort.bed
    sort -T $SORT_TMP -k1,1 -k2,2n -k3,3n $OUT/$sample/$sample.reads.undedup.summary_score.unsort.bed > $OUT/$sample/$sample.reads.undedup.summary_score.bed

    bed_inf=$OUT/$sample/${sample}.reads.undedup.summary_score.bed
    bed_outf=$OUT/$sample/${sample}.reads.undedup.summary_score.deduplicated.bed
    echo $bed_inf $bed_outf|python $SCRIPT_DIR/tools/deduplication_readsStartEnd.py
}

function merge_info_and_plot {
    IN=$HM/02-merge-outputs/info
    OUT=$HM/02-merge-outputs/info_merge
    
    mkdir -p $OUT/$sample

    echo $IN/$sample $OUT/$sample $sample.call.o.txt | python $SCRIPT_DIR/tools/merge_info.py

    # 注意：plot.R 在最后统一执行，这里不需要调用

}

function call_ratio {
    IN=$HM/02-merge-outputs/bedfiles
    OUT=$HM/02-merge-outputs/ratio
    chrom_sizes=${chromsizes[$sp]}

    reads=$IN/$sample/$sample.reads.undedup.summary_score.deduplicated.bed
    mkdir -p $OUT/$sample
    
    cat $reads |python $SCRIPT_DIR/tools/read_to_converted_site.py $chrom_sizes | grep -v NA  > $OUT/$sample/${sample}_sites.txt
    
    echo -e "Chromosome\tStart\tEnd\tFeature\tratio" > $OUT/$sample/${sample}_ratio.igv
    echo -e "Chromosome\tStart\tEnd\tFeature\tdepth" > $OUT/$sample/${sample}_depth.igv
    
    cat $OUT/$sample/${sample}_sites.txt | awk 'BEGIN{OFS="\t"}{if($4+$5>5) print $1,$2,$3,".",$6}' >> $OUT/${sample}/${sample}_ratio.igv
    cat $OUT/$sample/${sample}_sites.txt | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,".",$4+$5}' >> $OUT/${sample}/${sample}_depth.igv
    
    igv_sp=${igv_tdf_array[$sp]}
    
    igvtools toTDF -z 10 $OUT/${sample}/${sample}_depth.igv $OUT/${sample}/${sample}_depth.igv.tdf $igv_sp
    igvtools toTDF -z 10 $OUT/${sample}/${sample}_ratio.igv $OUT/${sample}/${sample}_ratio.igv.tdf $igv_sp
    
}


merge_bam_without_dedup

call_peak

# merge_bed_and_dedup

# 这次使用首尾去重
merge_bed_and_dedup_readsStartEnd

merge_info_and_plot

call_ratio

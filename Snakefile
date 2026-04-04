#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Bismark Pipeline with Snakemake
自动化bismark测序数据处理workflow
Author: Based on scripts by Shen Ke
Date: 2025
"""

import os
from pathlib import Path
import glob

# ============= 加载配置 =============
configfile: "config.yaml"

# ============= 全局变量 =============
# 获取Snakefile所在的绝对路径作为工作目录
WORK_DIR = os.path.dirname(os.path.abspath(workflow.snakefile))

HM = os.path.abspath(config["base_dir"])
SPECIES = config["species"]
N_PARTS = config["n_parts"]
SCRIPTS_DIR = config["paths"]["scripts"]
BISMARK_ENV = config["paths"]["bismark_env"]
RAW_DATA_SOURCE = os.path.abspath(config["raw_data_source"])
SAMPLE_RENAME_SCRIPT = os.path.abspath(config["paths"]["sample_rename_script"])

# ============= 自动检测样本 =============
RAWDATA_RAW_DIR = f"{HM}/00-rawdata_raw"
SAMPLES = []

# 优先从原始数据源检测样本（通过 _R1 文件推断样本名）
if os.path.exists(RAW_DATA_SOURCE):
    r1_files = glob.glob(f"{RAW_DATA_SOURCE}/*_R1*.fastq.gz") + \
               glob.glob(f"{RAW_DATA_SOURCE}/*_R1*.fq.gz") + \
               glob.glob(f"{RAW_DATA_SOURCE}/*.R1*.fastq.gz") + \
               glob.glob(f"{RAW_DATA_SOURCE}/*.R1*.fq.gz")
    
    for r1_file in r1_files:
        basename = os.path.basename(r1_file)
        # 提取样本名（去除 _R1 或 .R1 及后续部分）
        sample = basename.split('_R1')[0].split('.R1')[0]
        if sample and sample not in SAMPLES:
            SAMPLES.append(sample)
    
    if SAMPLES:
        print(f"[Snakemake] Detected {len(SAMPLES)} samples from source: {SAMPLES}")
    else:
        print(f"[Warning] No R1 fastq.gz files found in: {RAW_DATA_SOURCE}")
elif os.path.exists(RAWDATA_RAW_DIR):
    # 如果原始数据源不存在，从 00-rawdata_raw 检测（用于已运行过的情况）
    SAMPLES = [d for d in os.listdir(RAWDATA_RAW_DIR) 
               if os.path.isdir(os.path.join(RAWDATA_RAW_DIR, d)) and not d.startswith('.')]
    print(f"[Snakemake] Detected {len(SAMPLES)} samples from rawdata_raw: {SAMPLES}")
else:
    print(f"[Error] Neither raw_data_source ({RAW_DATA_SOURCE}) nor {RAWDATA_RAW_DIR} exists!")
    print("[Error] Pipeline cannot determine samples. Please check config.yaml")

# 样本检测完成
_ = None  # 结束 Python 代码块

# ============= 通配符约束 =============
wildcard_constraints:
    sample="[^/]+",
    part=r"\d+"

# ============= 辅助函数 =============
def get_all_part_bams(wildcards):
    """获取某个样本的所有part的bam文件"""
    return expand(f"{HM}/01-outputs/align/{{sample}}-part{{part}}/{{sample}}-part{{part}}.sorted.bam",
                  sample=wildcards.sample, part=range(N_PARTS))

def get_all_part_beds(wildcards):
    """获取某个样本的所有part的bed文件"""
    return expand(f"{HM}/02-merge-outputs/split_bedfiles/{{sample}}/{{sample}}-part{{part}}.reads.undedup.summary_score.bed",
                  sample=wildcards.sample, part=range(N_PARTS))

def get_all_part_infos(wildcards):
    """获取某个样本的所有part的info文件"""
    return expand(f"{HM}/02-merge-outputs/info/{{sample}}/{{sample}}-part{{part}}.call.o.txt",
                  sample=wildcards.sample, part=range(N_PARTS))

# ============= 最终输出目标 =============
def get_all_targets():
    """根据是否有样本返回不同的目标文件列表"""
    if not SAMPLES:
        # 如果没有样本，至少需要运行setup_rawdata来检查
        print("[Warning] No samples detected - will only run setup step")
        return [f"{HM}/00-rawdata_raw/.setup_done"]
    
    targets = []
    # 每个样本的最终sites文件
    targets.extend(expand(f"{HM}/02-merge-outputs/ratio/{{sample}}/{{sample}}_sites.txt", 
                          sample=SAMPLES))
    # 每个样本的IGV可视化文件
    targets.extend(expand(f"{HM}/02-merge-outputs/ratio/{{sample}}/{{sample}}_ratio.igv.tdf", 
                          sample=SAMPLES))
    targets.extend(expand(f"{HM}/02-merge-outputs/ratio/{{sample}}/{{sample}}_depth.igv.tdf", 
                          sample=SAMPLES))
    # 转化率图
    targets.append(f"{HM}/converted_ratio.pdf")
    
    return targets

rule all:
    input:
        get_all_targets()

# ============= 规则定义 =============

# 规则0: 初始化 - 将原始fastq文件链接到00-rawdata_raw目录
rule setup_rawdata:
    output:
        flag=f"{HM}/00-rawdata_raw/.setup_done",
        # 明确列出所有样本的R1和R2文件作为输出
        r1_files=expand(f"{HM}/00-rawdata_raw/{{sample}}/{{sample}}_R1.fastq.gz",
                       sample=SAMPLES),
        r2_files=expand(f"{HM}/00-rawdata_raw/{{sample}}/{{sample}}_R2.fastq.gz",
                       sample=SAMPLES)
    params:
        raw_source=RAW_DATA_SOURCE,
        out_dir=f"{HM}/00-rawdata_raw",
        script=SAMPLE_RENAME_SCRIPT,
        samples=" ".join(SAMPLES)  # 将SAMPLES列表转换为空格分隔的字符串
    log:
        f"{HM}/logs/setup/setup_rawdata.log"
    threads: config["resources"]["setup_rawdata"]["threads"]
    resources:
        mem_gb=config["resources"]["setup_rawdata"]["mem_gb"],
        partition=config["resources"]["setup_rawdata"]["partition"],
        time=config["resources"]["setup_rawdata"]["time"]
    shell:
        """
        # 创建日志目录和输出目录
        mkdir -p $(dirname {log})
        mkdir -p {params.out_dir}
        
        echo "[$(date)] Starting setup_rawdata" > {log}
        echo "Raw data source: {params.raw_source}" >> {log}
        echo "Output directory: {params.out_dir}" >> {log}
        echo "Samples to process: {params.samples}" >> {log}
        
        # 执行Sample_rename.new.sh脚本
        IN="{params.raw_source}" OUT="{params.out_dir}" bash {params.script} &>> {log}
        
        # 验证所有文件是否创建成功
        for sample in {params.samples}; do
            if [ ! -f "{params.out_dir}/$sample/${{sample}}_R1.fastq.gz" ]; then
                echo "Error: Failed to create R1 file for sample $sample" >> {log}
                exit 1
            fi
            if [ ! -f "{params.out_dir}/$sample/${{sample}}_R2.fastq.gz" ]; then
                echo "Error: Failed to create R2 file for sample $sample" >> {log}
                exit 1
            fi
        done
        
        # 创建完成标记
        touch {output.flag}
        
        echo "[$(date)] Completed setup_rawdata" >> {log}
        """

# 规则1: 拆分fastq文件
rule split_fastq:
    input:
        setup_flag=f"{HM}/00-rawdata_raw/.setup_done",
        r1=f"{HM}/00-rawdata_raw/{{sample}}/{{sample}}_R1.fastq.gz",
        r2=f"{HM}/00-rawdata_raw/{{sample}}/{{sample}}_R2.fastq.gz"
    output:
        # 使用标记文件表示拆分完成
        flag=f"{HM}/00-rawdata/.split_done/{{sample}}.done",
        # 明确列出所有part文件作为输出
        parts_r1=expand(f"{HM}/00-rawdata/{{{{sample}}}}-part{{part}}/{{{{sample}}}}-part{{part}}_R1.fastq.gz",
                        part=range(N_PARTS)),
        parts_r2=expand(f"{HM}/00-rawdata/{{{{sample}}}}-part{{part}}/{{{{sample}}}}-part{{part}}_R2.fastq.gz",
                        part=range(N_PARTS))
    params:
        outdir=f"{HM}/00-rawdata",
        sample="{sample}",
        n_parts=N_PARTS,
        input_dir=f"{HM}/00-rawdata_raw/{{sample}}",
        split_tool=f"{WORK_DIR}/tools/split_string.py"
    log:
        f"{HM}/logs/split/{{sample}}.log"
    threads: config["resources"]["split_fastq"]["threads"]
    resources:
        mem_gb=config["resources"]["split_fastq"]["mem_gb"],
        partition=config["resources"]["split_fastq"]["partition"],
        time=config["resources"]["split_fastq"]["time"]
    shell:
        """
        export PATH={BISMARK_ENV}:$PATH
        
        # 创建日志目录
        mkdir -p $(dirname {log})
        
        echo "[$(date)] Starting split for sample {wildcards.sample}" > {log}
        
        # 创建输出目录
        mkdir -p {params.outdir}
        mkdir -p $(dirname {output.flag})
        
        # 创建part目录
        for i in $(seq 0 $(({params.n_parts}-1))); do
            mkdir -p {params.outdir}/{params.sample}-part$i
        done
        
        # 生成并执行拆分命令 - R1
        split_cmd_r1=$(python {params.split_tool} -N {params.n_parts} \
                       -n {params.sample} -R R1 \
                       -i {params.input_dir} \
                       -o {params.outdir})
        echo "[$(date)] Splitting R1..." >> {log}
        echo "Command: $split_cmd_r1" >> {log}
        eval $split_cmd_r1 2>> {log}
        
        # 生成并执行拆分命令 - R2
        split_cmd_r2=$(python {params.split_tool} -N {params.n_parts} \
                       -n {params.sample} -R R2 \
                       -i {params.input_dir} \
                       -o {params.outdir})
        echo "[$(date)] Splitting R2..." >> {log}
        echo "Command: $split_cmd_r2" >> {log}
        eval $split_cmd_r2 2>> {log}
        
        # 创建完成标记
        touch {output.flag}
        echo "[$(date)] Split completed for {wildcards.sample}" >> {log}
        """

# 规则2: 处理每个part（质控、比对、提取突变）
rule process_part:
    input:
        # 依赖拆分完成标记
        flag=f"{HM}/00-rawdata/.split_done/{{sample}}.done",
        # 实际输入文件
        r1=f"{HM}/00-rawdata/{{sample}}-part{{part}}/{{sample}}-part{{part}}_R1.fastq.gz",
        r2=f"{HM}/00-rawdata/{{sample}}-part{{part}}/{{sample}}-part{{part}}_R2.fastq.gz"
    output:
        # 比对输出
        bam=f"{HM}/01-outputs/align/{{sample}}-part{{part}}/{{sample}}-part{{part}}.sorted.bam",
        bai=f"{HM}/01-outputs/align/{{sample}}-part{{part}}/{{sample}}-part{{part}}.sorted.bam.bai",
        # bed输出
        bed=f"{HM}/02-merge-outputs/split_bedfiles/{{sample}}/{{sample}}-part{{part}}.reads.undedup.summary_score.bed",
        # 信息输出
        info=f"{HM}/02-merge-outputs/info/{{sample}}/{{sample}}-part{{part}}.call.o.txt",
        # unmapped reads输出
        unmapped_r1=f"{HM}/01-outputs/unmapped/{{sample}}-part{{part}}/{{sample}}-part{{part}}_unmapped_R1.fq.gz",
        unmapped_r2=f"{HM}/01-outputs/unmapped/{{sample}}-part{{part}}/{{sample}}-part{{part}}_unmapped_R2.fq.gz"
    params:
        sample="{sample}",
        sample_part="{sample}-part{part}",
        species=SPECIES,
        hm=HM,
        process_script=f"{WORK_DIR}/do_process.sh",
        work_dir=WORK_DIR
    log:
        f"{HM}/logs/process/{{sample}}/{{sample}}-part{{part}}.log"
    threads: config["resources"]["process_part"]["threads"]
    resources:
        mem_gb=config["resources"]["process_part"]["mem_gb"],
        partition=config["resources"]["process_part"]["partition"],
        time=config["resources"]["process_part"]["time"]
    shell:
        """
        # 创建日志目录
        mkdir -p $(dirname {log})
        
        echo "[$(date)] Starting process_part for {params.sample_part}" > {log}
        
        # 执行do_process.sh脚本
        bash {params.process_script} {params.hm} {params.species} {threads} \
             {params.sample} {params.sample_part} &>> {log}
        
        echo "[$(date)] Completed process_part for {params.sample_part}" >> {log}
        """

# 规则3: 合并bam、call peaks、去重、计算转化率
rule merge_and_analysis:
    input:
        # 依赖：该样本的所有part的输出
        bams=get_all_part_bams,
        beds=get_all_part_beds,
        infos=get_all_part_infos
    output:
        # 合并后的bam
        merged_bam=f"{HM}/02-merge-outputs/align/{{sample}}/{{sample}}.sort_name.bam",
        # peaks
        peaks=f"{HM}/02-merge-outputs/bedfiles/{{sample}}/{{sample}}/{{sample}}_peaks.narrowPeak",
        high_peaks=f"{HM}/02-merge-outputs/bedfiles/{{sample}}/{{sample}}/high_score_peaks.bed",
        coverage_bw=f"{HM}/02-merge-outputs/bedfiles/{{sample}}/{{sample}}/coverage.bw",
        # 去重后的bed
        dedup_bed=f"{HM}/02-merge-outputs/bedfiles/{{sample}}/{{sample}}.reads.undedup.summary_score.deduplicated.bed",
        # 合并的info
        merged_info=f"{HM}/02-merge-outputs/info_merge/{{sample}}/{{sample}}.call.o.txt",
        # 最终ratio文件
        sites=f"{HM}/02-merge-outputs/ratio/{{sample}}/{{sample}}_sites.txt",
        ratio_igv=f"{HM}/02-merge-outputs/ratio/{{sample}}/{{sample}}_ratio.igv",
        ratio_tdf=f"{HM}/02-merge-outputs/ratio/{{sample}}/{{sample}}_ratio.igv.tdf",
        depth_tdf=f"{HM}/02-merge-outputs/ratio/{{sample}}/{{sample}}_depth.igv.tdf"
    params:
        sample="{sample}",
        species=SPECIES,
        hm=HM,
        merge_script=f"{WORK_DIR}/mergeBam_callPeak_dedup_callRatio.sh"
    log:
        f"{HM}/logs/merge/{{sample}}.log"
    threads: config["resources"]["merge_and_analysis"]["threads"]
    resources:
        mem_gb=config["resources"]["merge_and_analysis"]["mem_gb"],
        partition=config["resources"]["merge_and_analysis"]["partition"],
        time=config["resources"]["merge_and_analysis"]["time"]
    shell:
        """
        # 创建日志目录
        mkdir -p $(dirname {log})
        
        echo "[$(date)] Starting merge_and_analysis for {params.sample}" > {log}
        
        # 执行mergeBam_callPeak_dedup_callRatio.sh脚本
        bash {params.merge_script} {params.hm} {params.species} \
             {threads} {params.sample} &>> {log}
        
        echo "[$(date)] Completed merge_and_analysis for {params.sample}" >> {log}
        """

# 规则4: 绘制转化率图
rule plot_conversion:
    input:
        # 依赖所有样本的merged info文件
        expand(f"{HM}/02-merge-outputs/info_merge/{{sample}}/{{sample}}.call.o.txt",
               sample=SAMPLES)
    output:
        f"{HM}/converted_ratio.pdf"
    params:
        plot_script=f"{WORK_DIR}/tools/plot.R"
    log:
        f"{HM}/logs/plot/conversion.log"
    threads: config["resources"]["plot_conversion"]["threads"]
    resources:
        mem_gb=config["resources"]["plot_conversion"]["mem_gb"],
        partition=config["resources"]["plot_conversion"]["partition"],
        time=config["resources"]["plot_conversion"]["time"]
    shell:
        """
        # 创建日志目录
        mkdir -p $(dirname {log})
        
        echo "[$(date)] Starting plot_conversion" > {log}
        
        # 切换到工作目录执行R脚本（plot.R可能依赖当前目录）
        cd {HM}/02-merge-outputs/info_merge && Rscript {params.plot_script} &>> {log}
        
        # 将生成的PDF移动到Snakemake声明的输出路径
        SRC_PDF="{HM}/02-merge-outputs/info_merge/converted_ratio.pdf"
        if [ ! -f "$SRC_PDF" ]; then
            echo "Error: converted_ratio.pdf not generated at $SRC_PDF" >> {log}
            exit 1
        fi
        mv -f "$SRC_PDF" {output}
        
        echo "[$(date)] Completed plot_conversion" >> {log}
        """

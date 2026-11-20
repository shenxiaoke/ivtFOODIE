#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Bismark Pipeline with Snakemake
Automated Bismark sequencing data processing workflow
Author: Based on scripts by Shen Ke
Date: 2025
"""

import os
from pathlib import Path
import glob

# ============= Load config =============
configfile: "config.yaml"

# ============= Global variables =============
# Working directory: absolute path of this Snakefile's directory
WORK_DIR = os.path.dirname(os.path.abspath(workflow.snakefile))

HM = os.path.abspath(config["base_dir"])
SPECIES = config["species"]
N_PARTS = config["n_parts"]
SCRIPTS_DIR = config["paths"]["scripts"]
BISMARK_ENV = config["paths"]["bismark_env"]
RAW_DATA_SOURCE = os.path.abspath(config["raw_data_source"])
SAMPLE_RENAME_SCRIPT = os.path.abspath(config["paths"]["sample_rename_script"])

# ============= Auto-detect samples =============
RAWDATA_RAW_DIR = f"{HM}/00-rawdata_raw"
SAMPLES = []

# Prefer detecting samples from raw data source (infer names from _R1 files)
if os.path.exists(RAW_DATA_SOURCE):
    r1_files = glob.glob(f"{RAW_DATA_SOURCE}/*_R1*.fastq.gz") + \
               glob.glob(f"{RAW_DATA_SOURCE}/*_R1*.fq.gz") + \
               glob.glob(f"{RAW_DATA_SOURCE}/*.R1*.fastq.gz") + \
               glob.glob(f"{RAW_DATA_SOURCE}/*.R1*.fq.gz")
    
    for r1_file in r1_files:
        basename = os.path.basename(r1_file)
        # Extract sample name (strip _R1 or .R1 and everything after)
        sample = basename.split('_R1')[0].split('.R1')[0]
        if sample and sample not in SAMPLES:
            SAMPLES.append(sample)
    
    if SAMPLES:
        print(f"[Snakemake] Detected {len(SAMPLES)} samples from source: {SAMPLES}")
    else:
        print(f"[Warning] No R1 fastq.gz files found in: {RAW_DATA_SOURCE}")
elif os.path.exists(RAWDATA_RAW_DIR):
    # If raw data source is missing, detect from 00-rawdata_raw (re-runs)
    SAMPLES = [d for d in os.listdir(RAWDATA_RAW_DIR) 
               if os.path.isdir(os.path.join(RAWDATA_RAW_DIR, d)) and not d.startswith('.')]
    print(f"[Snakemake] Detected {len(SAMPLES)} samples from rawdata_raw: {SAMPLES}")
else:
    print(f"[Error] Neither raw_data_source ({RAW_DATA_SOURCE}) nor {RAWDATA_RAW_DIR} exists!")
    print("[Error] Pipeline cannot determine samples. Please check config.yaml")

# Sample detection done
_ = None  # End Python code block

# ============= Wildcard constraints =============
wildcard_constraints:
    sample="[^/]+",
    part=r"\d+"

# ============= Helper functions =============
def get_all_part_bams(wildcards):
    """BAM files for all parts of a sample."""
    return expand(f"{HM}/01-outputs/align/{{sample}}-part{{part}}/{{sample}}-part{{part}}.sorted.bam",
                  sample=wildcards.sample, part=range(N_PARTS))

def get_all_part_beds(wildcards):
    """BED files for all parts of a sample."""
    return expand(f"{HM}/02-merge-outputs/split_bedfiles/{{sample}}/{{sample}}-part{{part}}.reads.undedup.summary_score.bed",
                  sample=wildcards.sample, part=range(N_PARTS))

def get_all_part_infos(wildcards):
    """Info files for all parts of a sample."""
    return expand(f"{HM}/02-merge-outputs/info/{{sample}}/{{sample}}-part{{part}}.call.o.txt",
                  sample=wildcards.sample, part=range(N_PARTS))

# ============= Final output targets =============
def get_all_targets():
    """Return target outputs depending on whether samples were detected."""
    if not SAMPLES:
        # With no samples, at least run setup_rawdata for checks
        print("[Warning] No samples detected - will only run setup step")
        return [f"{HM}/00-rawdata_raw/.setup_done"]
    
    targets = []
    # Per-sample final sites file
    targets.extend(expand(f"{HM}/02-merge-outputs/ratio/{{sample}}/{{sample}}_sites.txt", 
                          sample=SAMPLES))
    # Per-sample IGV visualization tracks
    targets.extend(expand(f"{HM}/02-merge-outputs/ratio/{{sample}}/{{sample}}_ratio.igv.tdf", 
                          sample=SAMPLES))
    targets.extend(expand(f"{HM}/02-merge-outputs/ratio/{{sample}}/{{sample}}_depth.igv.tdf", 
                          sample=SAMPLES))
    # Conversion ratio plot
    targets.append(f"{HM}/converted_ratio.pdf")
    
    return targets

rule all:
    input:
        get_all_targets()

# ============= Rule definitions =============

# Rule 0: init - symlink raw FASTQ into 00-rawdata_raw
rule setup_rawdata:
    output:
        flag=f"{HM}/00-rawdata_raw/.setup_done",
        # List all sample R1/R2 outputs explicitly
        r1_files=expand(f"{HM}/00-rawdata_raw/{{sample}}/{{sample}}_R1.fastq.gz",
                       sample=SAMPLES),
        r2_files=expand(f"{HM}/00-rawdata_raw/{{sample}}/{{sample}}_R2.fastq.gz",
                       sample=SAMPLES)
    params:
        raw_source=RAW_DATA_SOURCE,
        out_dir=f"{HM}/00-rawdata_raw",
        script=SAMPLE_RENAME_SCRIPT,
        samples=" ".join(SAMPLES)  # Join SAMPLES into a space-separated string
    log:
        f"{HM}/logs/setup/setup_rawdata.log"
    threads: config["resources"]["setup_rawdata"]["threads"]
    resources:
        mem_gb=config["resources"]["setup_rawdata"]["mem_gb"],
        partition=config["resources"]["setup_rawdata"]["partition"],
        time=config["resources"]["setup_rawdata"]["time"]
    shell:
        """
        # Create log and output directories
        mkdir -p $(dirname {log})
        mkdir -p {params.out_dir}
        
        echo "[$(date)] Starting setup_rawdata" > {log}
        echo "Raw data source: {params.raw_source}" >> {log}
        echo "Output directory: {params.out_dir}" >> {log}
        echo "Samples to process: {params.samples}" >> {log}
        
        # Run Sample_rename.new.sh
        IN="{params.raw_source}" OUT="{params.out_dir}" bash {params.script} &>> {log}
        
        # Verify all expected files were created
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
        
        # Create completion flag
        touch {output.flag}
        
        echo "[$(date)] Completed setup_rawdata" >> {log}
        """

# Rule 1: split FASTQ files
rule split_fastq:
    input:
        setup_flag=f"{HM}/00-rawdata_raw/.setup_done",
        r1=f"{HM}/00-rawdata_raw/{{sample}}/{{sample}}_R1.fastq.gz",
        r2=f"{HM}/00-rawdata_raw/{{sample}}/{{sample}}_R2.fastq.gz"
    output:
        # Flag file marks split completion
        flag=f"{HM}/00-rawdata/.split_done/{{sample}}.done",
        # List all part files as outputs explicitly
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
        
        # Create log directory
        mkdir -p $(dirname {log})
        
        echo "[$(date)] Starting split for sample {wildcards.sample}" > {log}
        
        # Create output directories
        mkdir -p {params.outdir}
        mkdir -p $(dirname {output.flag})
        
        # Create per-part directories
        for i in $(seq 0 $(({params.n_parts}-1))); do
            mkdir -p {params.outdir}/{params.sample}-part$i
        done
        
        # Build and run split command for R1
        split_cmd_r1=$(python {params.split_tool} -N {params.n_parts} \
                       -n {params.sample} -R R1 \
                       -i {params.input_dir} \
                       -o {params.outdir})
        echo "[$(date)] Splitting R1..." >> {log}
        echo "Command: $split_cmd_r1" >> {log}
        eval $split_cmd_r1 2>> {log}
        
        # Build and run split command for R2
        split_cmd_r2=$(python {params.split_tool} -N {params.n_parts} \
                       -n {params.sample} -R R2 \
                       -i {params.input_dir} \
                       -o {params.outdir})
        echo "[$(date)] Splitting R2..." >> {log}
        echo "Command: $split_cmd_r2" >> {log}
        eval $split_cmd_r2 2>> {log}
        
        # Create completion flag
        touch {output.flag}
        echo "[$(date)] Split completed for {wildcards.sample}" >> {log}
        """

# Rule 2: process each part (QC, align, call conversions)
rule process_part:
    input:
        # Depends on split completion flag
        flag=f"{HM}/00-rawdata/.split_done/{{sample}}.done",
        # Actual input FASTQs
        r1=f"{HM}/00-rawdata/{{sample}}-part{{part}}/{{sample}}-part{{part}}_R1.fastq.gz",
        r2=f"{HM}/00-rawdata/{{sample}}-part{{part}}/{{sample}}-part{{part}}_R2.fastq.gz"
    output:
        # Alignment output
        bam=f"{HM}/01-outputs/align/{{sample}}-part{{part}}/{{sample}}-part{{part}}.sorted.bam",
        bai=f"{HM}/01-outputs/align/{{sample}}-part{{part}}/{{sample}}-part{{part}}.sorted.bam.bai",
        # BED output
        bed=f"{HM}/02-merge-outputs/split_bedfiles/{{sample}}/{{sample}}-part{{part}}.reads.undedup.summary_score.bed",
        # Info output
        info=f"{HM}/02-merge-outputs/info/{{sample}}/{{sample}}-part{{part}}.call.o.txt",
        # Unmapped reads output
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
        # Create log directory
        mkdir -p $(dirname {log})
        
        echo "[$(date)] Starting process_part for {params.sample_part}" > {log}
        
        # Run do_process.sh
        bash {params.process_script} {params.hm} {params.species} {threads} \
             {params.sample} {params.sample_part} &>> {log}
        
        echo "[$(date)] Completed process_part for {params.sample_part}" >> {log}
        """

# Rule 3: merge BAMs, call peaks, deduplicate, compute conversion ratio
rule merge_and_analysis:
    input:
        # All per-part outputs for this sample
        bams=get_all_part_bams,
        beds=get_all_part_beds,
        infos=get_all_part_infos
    output:
        # Merged BAM
        merged_bam=f"{HM}/02-merge-outputs/align/{{sample}}/{{sample}}.sort_name.bam",
        # peaks
        peaks=f"{HM}/02-merge-outputs/bedfiles/{{sample}}/{{sample}}/{{sample}}_peaks.narrowPeak",
        high_peaks=f"{HM}/02-merge-outputs/bedfiles/{{sample}}/{{sample}}/high_score_peaks.bed",
        coverage_bw=f"{HM}/02-merge-outputs/bedfiles/{{sample}}/{{sample}}/coverage.bw",
        # Deduplicated BED
        dedup_bed=f"{HM}/02-merge-outputs/bedfiles/{{sample}}/{{sample}}.reads.undedup.summary_score.deduplicated.bed",
        # Merged info
        merged_info=f"{HM}/02-merge-outputs/info_merge/{{sample}}/{{sample}}.call.o.txt",
        # Final ratio-related files
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
        # Create log directory
        mkdir -p $(dirname {log})
        
        echo "[$(date)] Starting merge_and_analysis for {params.sample}" > {log}
        
        # Run mergeBam_callPeak_dedup_callRatio.sh
        bash {params.merge_script} {params.hm} {params.species} \
             {threads} {params.sample} &>> {log}
        
        echo "[$(date)] Completed merge_and_analysis for {params.sample}" >> {log}
        """

# Rule 4: plot conversion ratio
rule plot_conversion:
    input:
        # Merged info for all samples
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
        # Create log directory
        mkdir -p $(dirname {log})
        
        echo "[$(date)] Starting plot_conversion" > {log}
        
        # Run R script from info_merge (plot.R may depend on CWD)
        cd {HM}/02-merge-outputs/info_merge && Rscript {params.plot_script} &>> {log}
        
        # Move generated PDF to Snakemake-declared output path
        SRC_PDF="{HM}/02-merge-outputs/info_merge/converted_ratio.pdf"
        if [ ! -f "$SRC_PDF" ]; then
            echo "Error: converted_ratio.pdf not generated at $SRC_PDF" >> {log}
            exit 1
        fi
        mv -f "$SRC_PDF" {output}
        
        echo "[$(date)] Completed plot_conversion" >> {log}
        """

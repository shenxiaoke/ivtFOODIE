"""
Motif discovery and scanning rules

功能：
    1. STREME motif富集：从de novo footprints中富集motif
    2. FIMO motif扫描：在基因组上查找TF结合位点
    3. 开放区交集：使用bedtools获取开放区的结合位点
"""

# ==================================================
# 获取输出目录（在规则定义时计算）
# ==================================================

OUTPUT_DIR = get_output_dir()

# ==================================================
# 辅助函数：从MEME文件提取motif列表
# ==================================================

def get_motifs_from_meme(meme_file):
    """
    从MEME文件中提取所有motif名称
    
    参数：
        meme_file: MEME文件路径
    
    返回：
        motif名称列表
    """
    import sys
    from pathlib import Path
    
    # 获取workflow目录
    workflow_dir = Path(workflow.basedir) if hasattr(workflow, 'basedir') else Path(workflow.snakefile).parent
    extract_script = workflow_dir / "scripts" / "motif" / "extract_motif_list.py"
    
    # 执行脚本提取motif列表
    import subprocess
    result = subprocess.run(
        [sys.executable, str(extract_script), meme_file],
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        return []
    
    motifs = result.stdout.strip().split(',')
    return [m for m in motifs if m]  # 过滤空字符串

def get_input_path_footprint(wildcards):
    """获取footprint bed文件的输入路径"""
    return f"{get_output_dir()}/{wildcards.tf}/04_footprints/{wildcards.tf}.footprint.{wildcards.method}.{wildcards.kd_thres_type}.bed"

def get_output_path_streme_meme(wildcards):
    """获取STREME meme文件的输出路径"""
    return f"{get_output_dir()}/{wildcards.tf}/05_motifs/streme/{wildcards.tf}.footprint.{wildcards.method}.{wildcards.kd_thres_type}.meme"

def get_output_path_streme_dir(wildcards):
    """获取STREME输出目录的路径"""
    return directory(f"{get_output_dir()}/{wildcards.tf}/05_motifs/streme/{wildcards.tf}.footprint.{wildcards.method}.{wildcards.kd_thres_type}.streme_out")

def get_input_path_streme_meme(wildcards):
    """获取STREME meme文件的输入路径"""
    return f"{get_output_dir()}/{wildcards.tf}/05_motifs/streme/{wildcards.tf}.footprint.{wildcards.method}.{wildcards.kd_thres_type}.meme"

def get_output_path_fimo_bed(wildcards):
    """获取FIMO bed文件的输出路径"""
    return f"{get_output_dir()}/{wildcards.tf}/05_motifs/fimo/{wildcards.tf}.footprint.{wildcards.method}.{wildcards.kd_thres_type}.motif_{wildcards.motif_name}.bed"

def get_output_path_fimo_dir(wildcards):
    """获取FIMO输出目录的路径"""
    return directory(f"{get_output_dir()}/{wildcards.tf}/05_motifs/fimo/{wildcards.tf}.footprint.{wildcards.method}.{wildcards.kd_thres_type}.motif_{wildcards.motif_name}.fimo_out")

def get_input_path_fimo_bed(wildcards):
    """获取FIMO bed文件的输入路径"""
    return f"{get_output_dir()}/{wildcards.tf}/05_motifs/fimo/{wildcards.tf}.footprint.{wildcards.method}.{wildcards.kd_thres_type}.motif_{wildcards.motif_name}.bed"

def get_output_path_intersect(wildcards):
    """获取intersect bed文件的输出路径"""
    return f"{get_output_dir()}/{wildcards.tf}/06_binding_sites/{wildcards.tf}.footprint.{wildcards.method}.{wildcards.kd_thres_type}.motif_{wildcards.motif_name}.open_region.bed"

# ==================================================
# Rule: STREME motif enrichment (CHECKPOINT)
# ==================================================
# 使用checkpoint来处理动态输出（motif数量未知）

checkpoint streme_enrich:
    input:
        footprint=get_input_path_footprint
    output:
        meme=OUTPUT_DIR + "/{tf}/05_motifs/streme/{tf}.footprint.{method}.{kd_thres_type}.meme",
        streme_dir=directory(OUTPUT_DIR + "/{tf}/05_motifs/streme/{tf}.footprint.{method}.{kd_thres_type}.streme_out")
    conda:
        "../envs/meme.yaml"
    threads: 1
    resources:
        mem_mb=16000,
        partition=get_cluster_partition,
        time=get_cluster_time
    log:
        "logs/streme/{tf}_{method}_{kd_thres_type}.log"
    benchmark:
        "benchmarks/streme/{tf}_{method}_{kd_thres_type}.txt"
    params:
        genome=lambda wildcards: config["genomes"][
            config["tfs"][wildcards.tf]["species"]
        ],
        minw=lambda wildcards: config.get("streme_minw", 6),
        maxw=lambda wildcards: config.get("streme_maxw", 15),
        nmotifs=lambda wildcards: config.get("streme_nmotifs", 10),
        thresh=lambda wildcards: config.get("streme_thresh", 0.05),
        totallength=lambda wildcards: config.get("streme_totallength", 4000000)
    shell:
        """
        set +e  # 允许命令失败而不中断
        python {workflow.basedir}/scripts/motif/streme_enrich.py \
            --footprint_bed {input.footprint} \
            --genome_fasta {params.genome} \
            --output_meme {output.meme} \
            --streme_dir {output.streme_dir} \
            --minw {params.minw} \
            --maxw {params.maxw} \
            --nmotifs {params.nmotifs} \
            --thresh {params.thresh} \
            --totallength {params.totallength} \
            > {log} 2>&1
        
        if [ $? -ne 0 ]; then
            echo "ERROR: STREME enrichment failed. Creating empty output files..." >> {log}
            mkdir -p $(dirname {output.meme})
            touch {output.meme}
            mkdir -p {output.streme_dir}
            echo "Empty output files created to allow workflow to continue." >> {log}
        fi
        """


# ==================================================
# Rule: FIMO motif scanning (REMOVED - using fimo_scan_per_motif instead)
# ==================================================
# Note: This rule has been removed. Use fimo_scan_per_motif rule instead.

# ==================================================
# Rule: Intersect with open regions (REMOVED - using intersect_open_regions_per_motif instead)
# ==================================================
# Note: This rule has been removed. Use intersect_open_regions_per_motif rule instead.


# ==================================================
# Rule: FIMO motif scanning per motif
# ==================================================

rule fimo_scan_per_motif:
    input:
        meme=get_input_path_streme_meme
    output:
        bed=OUTPUT_DIR + "/{tf}/05_motifs/fimo/{tf}.footprint.{method}.{kd_thres_type}.motif_{motif_name}.bed",
        fimo_dir=directory(OUTPUT_DIR + "/{tf}/05_motifs/fimo/{tf}.footprint.{method}.{kd_thres_type}.motif_{motif_name}.fimo_out")
    conda:
        "../envs/meme.yaml"
    threads: 1
    resources:
        mem_mb=32000,
        partition=get_cluster_partition,
        time=get_cluster_time
    log:
        "logs/fimo/{tf}_{method}_{kd_thres_type}_motif_{motif_name}.log"
    benchmark:
        "benchmarks/fimo/{tf}_{method}_{kd_thres_type}_motif_{motif_name}.txt"
    params:
        genome=lambda wildcards: config["genomes"][
            config["tfs"][wildcards.tf]["species"]
        ],
        thresh=lambda wildcards: config.get("fimo_thresh", 0.0001),
        max_stored=lambda wildcards: config.get("fimo_max_stored", 100000000),
        motif_name=lambda wildcards: wildcards.motif_name
    shell:
        """
        set +e  # 允许命令失败而不中断
        
        # 预检查：验证motif是否存在于MEME文件中
        echo "Pre-checking if motif '{params.motif_name}' exists in {input.meme}..." >> {log}
        motifs_list=$(python {workflow.basedir}/scripts/motif/extract_motif_list.py {input.meme} 2>> {log})
        if [ $? -ne 0 ]; then
            echo "WARNING: Failed to read MEME file {input.meme}" >> {log}
        else
            # 检查motif是否在列表中（使用逗号分隔的格式）
            if echo "$motifs_list" | grep -qE "(^|,){params.motif_name}(,|$)"; then
                echo "Motif '{params.motif_name}' found in MEME file. Proceeding with FIMO scan." >> {log}
            else
                echo "ERROR: Motif '{params.motif_name}' not found in {input.meme}" >> {log}
                echo "Available motifs: $motifs_list" >> {log}
                echo "Skipping FIMO scan for non-existent motif. This indicates a mismatch in the workflow." >> {log}
                # 不创建空文件，让workflow知道这个motif不存在
                exit 1
            fi
        fi
        
        # 运行FIMO扫描
        python {workflow.basedir}/scripts/motif/fimo_scan.py \
            --meme_file {input.meme} \
            --genome_fasta {params.genome} \
            --output_bed {output.bed} \
            --fimo_dir {output.fimo_dir} \
            --thresh {params.thresh} \
            --max_stored {params.max_stored} \
            --motif_name {params.motif_name} \
            > {log} 2>&1
        
        fimo_exit_code=$?
        
        # 检查是否是motif不存在的错误
        if [ $fimo_exit_code -ne 0 ]; then
            if grep -q "Motif.*not found" {log}; then
                echo "ERROR: Motif validation failed. This indicates a workflow configuration issue." >> {log}
                echo "The motif '{params.motif_name}' does not exist in {input.meme}" >> {log}
                echo "Please check the aggregate_motif_outputs rule to ensure it extracts motifs from the correct MEME file." >> {log}
                # 不创建空文件，让workflow知道这个motif不存在
                exit 1
            else
                echo "ERROR: FIMO scan failed for other reasons. Creating empty output files..." >> {log}
                mkdir -p $(dirname {output.bed})
                touch {output.bed}
                mkdir -p {output.fimo_dir}
                echo "Empty output files created to allow workflow to continue." >> {log}
            fi
        fi
        """


# ==================================================
# Rule: Intersect with open regions per motif
# ==================================================

rule intersect_open_regions_per_motif:
    input:
        fimo_bed=get_input_path_fimo_bed,
        open_region=lambda wildcards: config["open_regions"][
            config["tfs"][wildcards.tf]["cell_line"]
        ]
    output:
        OUTPUT_DIR + "/{tf}/06_binding_sites/{tf}.footprint.{method}.{kd_thres_type}.motif_{motif_name}.open_region.bed"
    conda:
        "../envs/base.yaml"
    threads: 1
    resources:
        mem_mb=8000,
        partition=get_cluster_partition,
        time=get_cluster_time
    log:
        "logs/intersect_open/{tf}_{method}_{kd_thres_type}_motif_{motif_name}.log"
    benchmark:
        "benchmarks/intersect_open/{tf}_{method}_{kd_thres_type}_motif_{motif_name}.txt"
    shell:
        """
        set +e  # 允许命令失败而不中断
        python {workflow.basedir}/scripts/motif/intersect_open_regions.py \
            --fimo_bed {input.fimo_bed} \
            --open_region {input.open_region} \
            --output {output} \
            > {log} 2>&1
        
        if [ $? -ne 0 ]; then
            echo "ERROR: Bedtools intersect failed. Creating empty output file..." >> {log}
            mkdir -p $(dirname {output})
            touch {output}
            echo "Empty output file created to allow workflow to continue." >> {log}
        fi
        """


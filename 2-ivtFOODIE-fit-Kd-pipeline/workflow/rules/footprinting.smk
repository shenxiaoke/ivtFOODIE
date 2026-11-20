"""
De novo footprinting规则：基于Kd阈值自动发现TF结合足迹

功能：
    对每种拟合方法（basic和flank_corrected）和每种阈值组合，
    生成独立的footprint结果
"""

# ==================================================
# 获取输出目录（在规则定义时计算）
# ==================================================

OUTPUT_DIR = get_output_dir()

# ==================================================
# 辅助函数
# ==================================================

def get_input_path_fit_kd(wildcards):
    """获取Kd拟合文件的输入路径"""
    return f"{get_output_dir()}/{wildcards.tf}/03_fitting/{wildcards.tf}.fittingKd.{wildcards.method}.tsv"

def get_output_path_footprint_bed(wildcards):
    """获取footprint bed文件的输出路径"""
    return f"{get_output_dir()}/{wildcards.tf}/04_footprints/{wildcards.tf}.footprint.{wildcards.method}.{wildcards.kd_thres_type}.bed"

def get_output_path_footprint_stats(wildcards):
    """获取footprint统计文件的输出路径"""
    return f"{get_output_dir()}/{wildcards.tf}/04_footprints/{wildcards.tf}.footprint.{wildcards.method}.{wildcards.kd_thres_type}.stats.txt"

# ==================================================
# Rule: De novo footprinting
# ==================================================

rule footprint_denovo:
    input:
        get_input_path_fit_kd
    output:
        bed=OUTPUT_DIR + "/{tf}/04_footprints/{tf}.footprint.{method}.{kd_thres_type}.bed",
        stats=OUTPUT_DIR + "/{tf}/04_footprints/{tf}.footprint.{method}.{kd_thres_type}.stats.txt"
    conda:
        "../envs/base.yaml"
    threads: 16  # 使用多进程并行
    resources:
        mem_mb=32000,
        partition=get_cluster_partition,
        time=get_cluster_time
    log:
        "logs/footprint_denovo/{tf}_{method}_{kd_thres_type}.log"
    benchmark:
        "benchmarks/footprint_denovo/{tf}_{method}_{kd_thres_type}.txt"
    params:
        seed_quantile=lambda wildcards: config["footprint_thresholds"][wildcards.kd_thres_type]["seed"],
        extend_quantile=lambda wildcards: config["footprint_thresholds"][wildcards.kd_thres_type]["extend"],
        min_len=lambda wildcards: config.get("footprint_min_len", 10),
        max_len=lambda wildcards: config.get("footprint_max_len", 30),
        max_gap=lambda wildcards: config.get("footprint_max_gap", 5),
        tolerance=lambda wildcards: config.get("footprint_tolerance", 0.1)
    shell:
        """
        {{
            # 运行footprinting
            python {DENOVO_FOOTPRINT} \
                --in_f {input} \
                --out_f {output.bed} \
                --seed_quantile {params.seed_quantile} \
                --extend_quantile {params.extend_quantile} \
                --num_processes {threads} \
                --min_len {params.min_len} \
                --max_len {params.max_len} \
                --max_gap {params.max_gap} \
                --tolerance {params.tolerance}

            # 生成统计信息
            echo "=== Footprinting Statistics ===" > {output.stats}
            echo "TF: {wildcards.tf}" >> {output.stats}
            echo "Method: {wildcards.method}" >> {output.stats}
            echo "Threshold type: {wildcards.kd_thres_type}" >> {output.stats}
            echo "Seed quantile: {params.seed_quantile}" >> {output.stats}
            echo "Extend quantile: {params.extend_quantile}" >> {output.stats}
            echo "Total footprints: $(wc -l < {output.bed})" >> {output.stats}
            echo "Date: $(date)" >> {output.stats}

        }} > {log} 2>&1
        """

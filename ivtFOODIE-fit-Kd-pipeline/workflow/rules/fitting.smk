"""
Kd拟合规则：从conversion ratio拟合TF-DNA结合的解离常数

功能：
    1. fit_kd_basic: 普通Kd拟合（无flank校正）
    2. fit_kd_flank: Flank校正Kd拟合（需要motif文件）
"""

# ==================================================
# 获取输出目录（在规则定义时计算）
# ==================================================

OUTPUT_DIR = get_output_dir()

# ==================================================
# 辅助函数
# ==================================================

def get_input_path_integrate(wildcards):
    """获取整合样本的输入路径"""
    return f"{get_output_dir()}/{wildcards.tf}/02_integrated/{wildcards.tf}.samples_ratio.tsv"

def get_output_path_fit_kd_basic(wildcards):
    """获取basic Kd拟合的输出路径"""
    return f"{get_output_dir()}/{wildcards.tf}/03_fitting/{wildcards.tf}.fittingKd.basic.tsv"

def get_output_path_fit_kd_flank(wildcards):
    """获取flank Kd拟合的输出路径"""
    return f"{get_output_dir()}/{wildcards.tf}/03_fitting/{wildcards.tf}.fittingKd.flank_corrected.tsv"

# ==================================================
# Rule: 普通Kd拟合（无flank校正）
# ==================================================

rule fit_kd_basic:
    input:
        get_input_path_integrate
    output:
        OUTPUT_DIR + "/{tf}/03_fitting/{tf}.fittingKd.basic.tsv"
    conda:
        "../envs/base.yaml"
    threads: 1
    resources:
        mem_mb=16000,
        partition=get_cluster_partition,
        time=get_cluster_time
    log:
        "logs/fit_kd_basic/{tf}.log"
    benchmark:
        "benchmarks/fit_kd_basic/{tf}.txt"
    params:
        concentrations=lambda wildcards: ",".join(
            map(str, get_sample_info(wildcards.tf)["concentrations"])
        ),
        fitting_indices=lambda wildcards: ",".join(
            map(str, get_sample_info(wildcards.tf)["fitting_indices"])
        ),
        fix_n_to_1=lambda wildcards: str(config["tfs"][wildcards.tf].get("fix_n_to_1", True)),
        pc_name="ratio_adjust"  # 固定使用ratio_adjust（由integration脚本生成）
    shell:
        """
        python {KD_FITTER} \
            --in_f {input} \
            --out_f {output} \
            --concentration {params.concentrations} \
            --fitting_used_index {params.fitting_indices} \
            --fix_n_to_1 {params.fix_n_to_1} \
            --PC_name {params.pc_name} \
            > {log} 2>&1
        """


# ==================================================
# Rule: Flank校正Kd拟合（需要preliminary motif）
# ==================================================

rule fit_kd_flank:
    input:
        ratio=get_input_path_integrate,
        # preliminary motif: 从basic footprint结果中提取
        # 或者用户提供的已知motif文件
        motif=lambda wildcards: (
            config["tfs"][wildcards.tf].get("preliminary_motif") or
            get_preliminary_motif(wildcards.tf)
        )
    output:
        OUTPUT_DIR + "/{tf}/03_fitting/{tf}.fittingKd.flank_corrected.tsv"
    conda:
        "../envs/base.yaml"
    threads: 1
    resources:
        mem_mb=16000,
        partition=get_cluster_partition,
        time=get_cluster_time
    log:
        "logs/fit_kd_flank/{tf}.log"
    benchmark:
        "benchmarks/fit_kd_flank/{tf}.txt"
    params:
        concentrations=lambda wildcards: ",".join(
            map(str, get_sample_info(wildcards.tf)["concentrations"])
        ),
        fitting_indices=lambda wildcards: ",".join(
            map(str, get_sample_info(wildcards.tf)["fitting_indices"])
        ),
        fix_n_to_1=lambda wildcards: str(config["tfs"][wildcards.tf].get("fix_n_to_1", True)),
        pc_name="ratio_adjust",  # 固定使用ratio_adjust（由integration脚本生成）
        flank_size=lambda wildcards: config.get("flank_size", 50),
        motif_span=lambda wildcards: config.get("motif_span", 100),
        motif_pick_ratio=lambda wildcards: config.get("motif_pick_ratio", 0.1)
    shell:
        """
        python {KD_FITTER} \
            --in_f {input.ratio} \
            --out_f {output} \
            --concentration {params.concentrations} \
            --fitting_used_index {params.fitting_indices} \
            --fix_n_to_1 {params.fix_n_to_1} \
            --motif_f {input.motif} \
            --PC_name {params.pc_name} \
            --motif_pick_ratio {params.motif_pick_ratio} \
            --flank_size {params.flank_size} \
            --motif_span {params.motif_span} \
            > {log} 2>&1
        """

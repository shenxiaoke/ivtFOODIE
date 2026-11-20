"""
样本整合规则：整合多个样本的ratio数据

功能：
    从多个样本的tabix压缩文件中提取数据，整合为单一DataFrame
"""

# ==================================================
# 获取输出目录（在规则定义时计算）
# ==================================================

OUTPUT_DIR = get_output_dir()

# ==================================================
# 辅助函数
# ==================================================

def get_output_path_integrate(wildcards):
    """获取整合样本的输出路径"""
    return f"{get_output_dir()}/{wildcards.tf}/02_integrated/{wildcards.tf}.samples_ratio.tsv"

def get_input_path_pc_control(wildcards):
    """获取PC对照文件的输入路径"""
    return f"{get_output_dir()}/{wildcards.tf}/01_compressed/PC_control_sites.txt.gz"

def get_input_path_pc_index(wildcards):
    """获取PC对照索引文件的输入路径"""
    return f"{get_output_dir()}/{wildcards.tf}/01_compressed/PC_control_sites.txt.gz.tbi"

# ==================================================
# Rule: 整合多样本ratio数据
# ==================================================

rule integrate_samples:
    input:
        conversion_qc=conversion_qc_gate,
        # 样本文件：动态根据config生成
        samples=lambda wildcards: expand(
            f"{get_output_dir()}/{{tf}}/01_compressed/{{sample}}_sites.txt.gz",
            tf=wildcards.tf,
            sample=get_sample_names(wildcards.tf)
        ),
        # 样本索引文件
        indices=lambda wildcards: expand(
            f"{get_output_dir()}/{{tf}}/01_compressed/{{sample}}_sites.txt.gz.tbi",
            tf=wildcards.tf,
            sample=get_sample_names(wildcards.tf)
        ),
        # PC对照文件（workflow自动压缩）
        pc_control=get_input_path_pc_control,
        pc_index=get_input_path_pc_index,
        # Open region BED文件
        open_region=lambda wildcards: config["open_regions"][
            config["tfs"][wildcards.tf]["cell_line"]
        ]
    output:
        OUTPUT_DIR + "/{tf}/02_integrated/{tf}.samples_ratio.tsv"
    conda:
        "../envs/base.yaml"
    threads: 1
    resources:
        mem_mb=16000,
        partition=get_cluster_partition,
        time=get_cluster_time
    log:
        "logs/integrate_samples/{tf}.log"
    benchmark:
        "benchmarks/integrate_samples/{tf}.txt"
    params:
        tf_name=lambda wildcards: wildcards.tf,
        sample_path=lambda wildcards: Path(get_output_dir()) / wildcards.tf / "01_compressed",
        min_depth=lambda wildcards: config.get("min_depth", 10),
        genome=lambda wildcards: config["genomes"][
            config["tfs"][wildcards.tf]["species"]
        ]
    script:
        "../scripts/integration/integrate_samples.py"

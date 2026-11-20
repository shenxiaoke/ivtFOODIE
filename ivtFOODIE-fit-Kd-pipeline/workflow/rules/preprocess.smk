"""
预处理规则：压缩和索引ratio文件

功能：
    1. 使用bgzip压缩.txt文件为.txt.gz
    2. 对已使用bgzip压缩并建立索引的.txt.gz文件，直接复用原始.gz和.tbi文件
    3. 使用tabix建立索引文件.txt.gz.tbi

注意：
    - tabix索引要求输入文件已按坐标排序
    - tabix参数说明：
      -s: 染色体列号（1-based）
      -b: 起始位置列号（1-based）
      -e: 结束位置列号（1-based，如果是单碱基则与-b相同）
"""

# ==================================================
# 规则优先级：解决PC对照和样本文件的歧义
# ==================================================

# PC对照文件优先使用专门的规则
ruleorder: compress_pc > compress_ratio
ruleorder: index_pc > index_ratio

# ==================================================
# 获取输出目录（在规则定义时计算）
# ==================================================

OUTPUT_DIR = get_output_dir()

# ==================================================
# 辅助函数：获取样本的原始txt或txt.gz文件路径
# ==================================================

def get_output_path_compress_ratio(wildcards):
    """获取压缩ratio文件的输出路径"""
    return f"{get_output_dir()}/{wildcards.tf}/01_compressed/{wildcards.sample}_sites.txt.gz"

def get_output_path_compress_pc(wildcards):
    """获取压缩PC文件的输出路径"""
    return f"{get_output_dir()}/{wildcards.tf}/01_compressed/PC_control_sites.txt.gz"

def get_output_path_index_ratio(wildcards):
    """获取索引ratio文件的输出路径"""
    return f"{get_output_dir()}/{wildcards.tf}/01_compressed/{wildcards.sample}_sites.txt.gz.tbi"

def get_output_path_index_pc(wildcards):
    """获取索引PC文件的输出路径"""
    return f"{get_output_dir()}/{wildcards.tf}/01_compressed/PC_control_sites.txt.gz.tbi"

def get_output_path_validate(wildcards):
    """获取验证标志文件的输出路径"""
    return f"{get_output_dir()}/{wildcards.tf}/01_compressed/{wildcards.sample}.validated"

def get_input_path_compress_ratio(wildcards):
    """获取压缩ratio文件的输入路径"""
    return f"{get_output_dir()}/{wildcards.tf}/01_compressed/{wildcards.sample}_sites.txt.gz"

def get_input_path_compress_pc(wildcards):
    """获取压缩PC文件的输入路径"""
    return f"{get_output_dir()}/{wildcards.tf}/01_compressed/PC_control_sites.txt.gz"

def get_sample_txt_path(wildcards):
    """
    获取指定样本的原始txt或txt.gz文件路径

    使用Snakefile中的get_sample_info函数获取样本信息
    """
    from pathlib import Path

    # 获取该TF的所有样本信息
    info = get_sample_info(wildcards.tf)
    sample_ids = info['sample_ids']
    sample_files = info['sample_files']

    # 查找匹配的样本
    try:
        idx = sample_ids.index(wildcards.sample)
        return sample_files[idx]
    except ValueError:
        raise ValueError(
            f"Sample '{wildcards.sample}' not found in TF '{wildcards.tf}' configuration"
        )

def get_pc_txt_path(wildcards):
    """获取PC对照的原始txt或txt.gz文件路径"""
    return config["tfs"][wildcards.tf]["pc_control"]

def get_existing_tabix_index(file_path):
    """如果输入已经是bgzip压缩文件，返回必须随附的tabix索引路径"""
    if str(file_path).endswith(".gz"):
        return f"{file_path}.tbi"
    return []

def get_sample_existing_tabix_index(wildcards):
    """获取已压缩样本输入随附的tabix索引路径"""
    return get_existing_tabix_index(get_sample_txt_path(wildcards))

def get_pc_existing_tabix_index(wildcards):
    """获取已压缩PC输入随附的tabix索引路径"""
    return get_existing_tabix_index(get_pc_txt_path(wildcards))


# ==================================================
# Rule: 压缩样本ratio文件
# ==================================================

rule compress_ratio:
    input:
        get_sample_txt_path
    output:
        OUTPUT_DIR + "/{tf}/01_compressed/{sample}_sites.txt.gz"
    conda:
        "../envs/genomics.yaml"
    threads: 2
    resources:
        partition=get_cluster_partition,
        time=get_cluster_time
    log:
        "logs/compress_ratio/{tf}_{sample}.log"
    benchmark:
        "benchmarks/compress_ratio/{tf}_{sample}.txt"
    shell:
        """
        case "{input}" in
            *.gz)
                ln -sfn "$(realpath "{input}")" "{output}" 2> "{log}"
                ;;
            *)
                bgzip -c -@ {threads} "{input}" > "{output}" 2> "{log}"
                ;;
        esac
        """


# ==================================================
# Rule: 压缩PC对照ratio文件
# ==================================================

rule compress_pc:
    input:
        get_pc_txt_path
    output:
        OUTPUT_DIR + "/{tf}/01_compressed/PC_control_sites.txt.gz"
    conda:
        "../envs/genomics.yaml"
    threads: 2
    resources:
        partition=get_cluster_partition,
        time=get_cluster_time
    log:
        "logs/compress_pc/{tf}.log"
    benchmark:
        "benchmarks/compress_pc/{tf}.txt"
    shell:
        """
        case "{input}" in
            *.gz)
                ln -sfn "$(realpath "{input}")" "{output}" 2> "{log}"
                ;;
            *)
                bgzip -c -@ {threads} "{input}" > "{output}" 2> "{log}"
                ;;
        esac
        """


# ==================================================
# Rule: 建立tabix索引（样本文件）
# ==================================================

rule index_ratio:
    input:
        gz=get_input_path_compress_ratio,
        existing_tbi=get_sample_existing_tabix_index
    output:
        OUTPUT_DIR + "/{tf}/01_compressed/{sample}_sites.txt.gz.tbi"
    conda:
        "../envs/genomics.yaml"
    resources:
        partition=get_cluster_partition,
        time=get_cluster_time
    log:
        "logs/index_ratio/{tf}_{sample}.log"
    benchmark:
        "benchmarks/index_ratio/{tf}_{sample}.txt"
    params:
        # ratio文件格式：chr, pos, ..., conversion, unconversion, ratio
        # 假设列顺序：1=chr, 2=pos
        sequence_col=1,   # 染色体列（1-based）
        start_col=2,      # 起始位置列（1-based）
        # 单碱基位点，所以end_col = start_col
    shell:
        """
        if [ -n "{input.existing_tbi}" ]; then
            ln -sfn "$(realpath "{input.existing_tbi}")" "{output}" 2> "{log}"
        else
            tabix -s {params.sequence_col} -b {params.start_col} -e {params.start_col} "{input.gz}" 2> "{log}"
        fi
        """


# ==================================================
# Rule: 建立tabix索引（PC对照文件）
# ==================================================

rule index_pc:
    input:
        gz=get_input_path_compress_pc,
        existing_tbi=get_pc_existing_tabix_index
    output:
        OUTPUT_DIR + "/{tf}/01_compressed/PC_control_sites.txt.gz.tbi"
    conda:
        "../envs/genomics.yaml"
    resources:
        partition=get_cluster_partition,
        time=get_cluster_time
    log:
        "logs/index_pc/{tf}.log"
    benchmark:
        "benchmarks/index_pc/{tf}.txt"
    params:
        sequence_col=1,
        start_col=2,
    shell:
        """
        if [ -n "{input.existing_tbi}" ]; then
            ln -sfn "$(realpath "{input.existing_tbi}")" "{output}" 2> "{log}"
        else
            tabix -s {params.sequence_col} -b {params.start_col} -e {params.start_col} "{input.gz}" 2> "{log}"
        fi
        """


# ==================================================
# Rule: 验证压缩和索引的完整性
# ==================================================

rule validate_compression:
    input:
        gz=get_input_path_compress_ratio,
        tbi=OUTPUT_DIR + "/{tf}/01_compressed/{sample}_sites.txt.gz.tbi"
    output:
        flag=OUTPUT_DIR + "/{tf}/01_compressed/{sample}.validated"
    conda:
        "../envs/genomics.yaml"
    resources:
        partition=get_cluster_partition,
        time=get_cluster_time
    log:
        "logs/validate_compression/{tf}_{sample}.log"
    shell:
        """
        {{
            echo "Validating compression and indexing for {wildcards.sample}..."

            # 验证gzip文件完整性
            gzip -t {input.gz}
            echo "gzip integrity: OK"

            # 验证tabix索引可用
            tabix -l {input.gz} > /dev/null
            echo "tabix index: OK"

            # 验证文件格式（至少有6列）
            zcat {input.gz} | head -1 | awk '{{if (NF < 6) exit 1}}'
            echo "file format: OK (at least 6 columns)"

            # 创建验证标志文件
            echo "Validation passed on $(date)" > {output.flag}

        }} > {log} 2>&1
        """

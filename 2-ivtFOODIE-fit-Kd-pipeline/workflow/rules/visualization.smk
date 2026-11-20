"""
Visualization rules

功能：
    转化率可视化：绘制不同浓度下的转化率图
"""

# ==================================================
# 获取输出目录（在规则定义时计算）
# ==================================================

OUTPUT_DIR = get_output_dir()

# ==================================================
# 辅助函数
# ==================================================

def get_input_path_intersect(wildcards):
    """获取intersect bed文件的输入路径"""
    return f"{get_output_dir()}/{wildcards.tf}/06_binding_sites/{wildcards.tf}.footprint.{wildcards.method}.{wildcards.kd_thres_type}.motif_{wildcards.motif_name}.open_region.bed"

def get_output_path_visualization(wildcards):
    """获取可视化文件的输出路径"""
    return f"{get_output_dir()}/{wildcards.tf}/07_visualization/{wildcards.tf}.footprint.{wildcards.method}.{wildcards.kd_thres_type}.motif_{wildcards.motif_name}.conversion_plot.pdf"

# ==================================================
# Rule: Motif-centered conversion rate plot (REMOVED - using motif_conversion_plot_per_motif instead)
# ==================================================
# Note: This rule has been removed. Use motif_conversion_plot_per_motif rule instead.


# ==================================================
# Rule: Motif-centered conversion rate plot per motif
# ==================================================

rule motif_conversion_plot_per_motif:
    input:
        binding_sites=get_input_path_intersect,
        samples=lambda wildcards: expand(
            f"{get_output_dir()}/{{tf}}/01_compressed/{{sample}}_sites.txt.gz",
            tf=wildcards.tf,
            sample=get_sample_names(wildcards.tf)
        ),
        pc_control=get_input_path_pc_control
    output:
        OUTPUT_DIR + "/{tf}/07_visualization/{tf}.footprint.{method}.{kd_thres_type}.motif_{motif_name}.conversion_plot.pdf"
    conda:
        "../envs/base.yaml"
    threads: 1
    resources:
        mem_mb=16000,
        partition=get_cluster_partition,
        time=get_cluster_time
    log:
        "logs/visualization/{tf}_{method}_{kd_thres_type}_motif_{motif_name}.log"
    benchmark:
        "benchmarks/visualization/{tf}_{method}_{kd_thres_type}_motif_{motif_name}.txt"
    params:
        tf_name=lambda wildcards: wildcards.tf,
        motif_pick_ratio=lambda wildcards: 1.0,  # 使用所有结合位点（因为已经是单个motif）
        min_depth=lambda wildcards: config.get("visualization", {}).get("min_depth_plot", 5),
        span=lambda wildcards: config.get("motif_span", 100),
        normalize_flank=lambda wildcards: config.get("visualization", {}).get("normalize_flank", True),
        plot_format=lambda wildcards: config.get("visualization", {}).get("plot_format", "pdf"),
        plot_dpi=lambda wildcards: config.get("visualization", {}).get("plot_dpi", 300),
        concentrations=lambda wildcards: ",".join(
            map(str, get_sample_info(wildcards.tf)["concentrations"])
        ),
        sample_ids=lambda wildcards: ",".join(get_sample_names(wildcards.tf))
    shell:
        """
        set +e  # 允许命令失败而不中断
        # 将 samples 列表转换为逗号分隔的字符串
        SAMPLES_STR=$(echo "{input.samples}" | sed 's/ /,/g')
        
        # 构建命令，根据 normalize_flank 参数决定是否添加标志
        NORMALIZE_FLAG=""
        if [ "{params.normalize_flank}" = "True" ]; then
            NORMALIZE_FLAG="--normalize_flank"
        fi
        
        python {workflow.basedir}/scripts/visualization/motif_conversion_plot.py \
            --binding_sites {input.binding_sites} \
            --samples "$SAMPLES_STR" \
            --sample_ids {params.sample_ids} \
            --concentrations {params.concentrations} \
            --pc_control {input.pc_control} \
            --output {output} \
            --motif_pick_ratio {params.motif_pick_ratio} \
            --min_depth {params.min_depth} \
            --span {params.span} \
            --plot_format {params.plot_format} \
            --plot_dpi {params.plot_dpi} \
            $NORMALIZE_FLAG \
            > {log} 2>&1
        
        if [ $? -ne 0 ]; then
            echo "ERROR: Motif conversion plot failed. Creating empty output file..." >> {log}
            mkdir -p $(dirname {output})
            touch {output}
            echo "Empty output file created to allow workflow to continue." >> {log}
        fi
        """


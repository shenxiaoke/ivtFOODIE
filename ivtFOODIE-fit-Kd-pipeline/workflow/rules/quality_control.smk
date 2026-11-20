"""Overview conversion plot is a prerequisite for sample integration."""
import json


def conversion_qc_enabled():
    return config.get("conversion_qc", {}).get("enabled", True)


def conversion_qc_dir():
    group = "-".join(sorted(get_all_tfs()))
    return str(Path(get_output_dir()) / "00_qc" / group)


def conversion_qc_gate(wildcards):
    if conversion_qc_enabled():
        return [conversion_qc_dir() + "/sampling_metadata.json"]
    return []


def conversion_qc_inputs(wildcards):
    files = []
    for tf in get_all_tfs():
        paths = [f"{get_output_dir()}/{tf}/01_compressed/PC_control_sites.txt.gz"]
        paths += [f"{get_output_dir()}/{tf}/01_compressed/{s}_sites.txt.gz" for s in get_sample_names(tf)]
        files.extend(paths)
        files.extend(path + ".tbi" for path in paths)
        files.append(config["open_regions"][config["tfs"][tf]["cell_line"]])
    return list(dict.fromkeys(files))


def conversion_qc_spec(wildcards):
    specification = {"options": config.get("conversion_qc", {}), "tfs": {}}
    for tf in get_all_tfs():
        info = get_sample_info(tf)
        specification["tfs"][tf] = {
            "pc_control": f"{get_output_dir()}/{tf}/01_compressed/PC_control_sites.txt.gz",
            "pc_source": config["tfs"][tf]["pc_control"],
            "open_regions": config["open_regions"][config["tfs"][tf]["cell_line"]],
            "samples": [
                {"sample_id": name, "concentration_M": float(concentration), "original_path": original,
                 "path": f"{get_output_dir()}/{tf}/01_compressed/{name}_sites.txt.gz"}
                for name, concentration, original in zip(info["sample_ids"], info["concentrations"], info["sample_files"])
            ],
        }
    return json.dumps(specification)


rule conversion_qc:
    input:
        conversion_qc_inputs
    output:
        png=conversion_qc_dir() + "/mean_conversion_vs_concentration.png",
        pdf=conversion_qc_dir() + "/mean_conversion_vs_concentration.pdf",
        summary=conversion_qc_dir() + "/mean_conversion.tsv",
        metadata=conversion_qc_dir() + "/sampling_metadata.json"
    conda:
        "../envs/base.yaml"
    threads: 1
    resources:
        mem_mb=lambda wc: config.get("conversion_qc", {}).get("mem_mb", 16000),
        partition=get_cluster_partition,
        time=get_cluster_time
    log:
        "logs/conversion_qc.log"
    params:
        spec=conversion_qc_spec,
        output_dir=conversion_qc_dir(),
        script=str(WORKFLOW_DIR / "scripts" / "quality_control" / "conversion_qc.py")
    shell:
        "python {params.script:q} --spec {params.spec:q} --output-dir {params.output_dir:q} > {log:q} 2>&1"


rule all_qc:
    input:
        conversion_qc_gate

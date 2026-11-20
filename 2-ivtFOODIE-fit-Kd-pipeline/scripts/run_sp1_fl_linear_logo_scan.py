#!/usr/bin/env python3

import argparse
import csv
import math
import shutil
import subprocess
import sys
from pathlib import Path

import yaml


ANALYSIS_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_ROOT = ANALYSIS_ROOT / "ivtFOODIE_postprocess_codeRefactoring"
RESULTS_DIR = ANALYSIS_ROOT / "results" / "SP1-FL"
CONFIG_PATH = WORKFLOW_ROOT / "config" / "config-SP1-FL.yaml"
LINEAR_LOGO_PATH = Path(
    "/dshare/xielab/analysisdata/TFC/shenk/analysis/WZ/BANC_FOODIE/260602_FOODIEin/"
    "seq2FL-DBD_linearKd_model/logos/unfiltered/SP1/flank_corrected/FL_neglog10Kd/ols/linear_logo.tsv"
)
PLOT_SCRIPT = (
    WORKFLOW_ROOT
    / "workflow"
    / "scripts"
    / "visualization"
    / "motif_conversion_plot.py"
)
PLOT_PYTHON = Path("/share/home/shenk/local/anaconda3/envs/bismark/bin/python")
BEDTOOLS = Path("/share/home/shenk/local/bedtools2.26.0/bedtools")
FIMO = Path("/share/home/shenk/meme/bin/fimo")


def parse_args():
    parser = argparse.ArgumentParser(
        description="One-off SP1-FL linear_logo motif scan and conversion plot."
    )
    parser.add_argument(
        "--start-position",
        type=int,
        default=10,
        help="Inclusive start position from linear_logo.tsv.",
    )
    parser.add_argument(
        "--end-position",
        type=int,
        default=20,
        help="Inclusive end position from linear_logo.tsv.",
    )
    parser.add_argument(
        "--fimo-thresh",
        type=float,
        default=1e-4,
        help="FIMO p-value threshold.",
    )
    parser.add_argument(
        "--output-tag",
        default="FL_neglog10Kd_ols_linear_logo_pos10_20",
        help="Output tag under results/SP1-FL.",
    )
    parser.add_argument(
        "--tf-name",
        default="SP1-FL",
        help="TF key in the config file, e.g. SP1 or SP1-FL.",
    )
    parser.add_argument(
        "--results-tf",
        default=None,
        help="Results subdirectory under results/. Defaults to --tf-name.",
    )
    parser.add_argument(
        "--config-path",
        type=Path,
        default=None,
        help="Config YAML. Defaults to config/config-<tf-name>.yaml.",
    )
    parser.add_argument(
        "--logo-path",
        type=Path,
        default=LINEAR_LOGO_PATH,
        help="linear_logo.tsv or coef_logo_zero_sum.tsv to convert to a MEME motif.",
    )
    parser.add_argument(
        "--motif-name",
        default=None,
        help="Motif name to write into MEME/BED outputs.",
    )
    parser.add_argument(
        "--binding-sites-bed",
        type=Path,
        default=None,
        help="Existing TFBS BED to plot directly; skips logo conversion and FIMO.",
    )
    parser.add_argument(
        "--top-n-tfbs",
        type=int,
        default=None,
        help="Use only the top N TFBS by BED score column for plotting.",
    )
    parser.add_argument(
        "--visualization-dir-name",
        default="07_visualization",
        help="Visualization subdirectory under the output tag.",
    )
    return parser.parse_args()


def read_config(path: Path):
    with path.open() as handle:
        return yaml.safe_load(handle)


def read_linear_logo(path: Path, start_pos: int, end_pos: int):
    rows_by_position = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            pos = int(row["position"])
            if pos < start_pos or pos > end_pos:
                continue
            rows_by_position.setdefault(pos, {})[row["base"]] = float(row["coef_logo_zero_sum"])

    expected_positions = list(range(start_pos, end_pos + 1))
    missing = [pos for pos in expected_positions if pos not in rows_by_position]
    if missing:
        raise ValueError(f"Missing positions in linear_logo.tsv: {missing}")

    for pos in expected_positions:
        missing_bases = [base for base in "ACGT" if base not in rows_by_position[pos]]
        if missing_bases:
            raise ValueError(f"Position {pos} missing bases: {missing_bases}")

    return rows_by_position


def pwm_from_logo(rows_by_position):
    pwm = []
    for pos in sorted(rows_by_position):
        weights = [10 ** rows_by_position[pos][base] for base in "ACGT"]
        total = sum(weights)
        probs = [weight / total for weight in weights]
        pwm.append((pos, probs))
    return pwm


def write_meme(meme_path: Path, motif_name: str, pwm):
    meme_path.parent.mkdir(parents=True, exist_ok=True)
    with meme_path.open("w") as handle:
        handle.write("MEME version 5\n\n")
        handle.write("ALPHABET= ACGT\n\n")
        handle.write("strands: + -\n\n")
        handle.write("Background letter frequencies\n")
        handle.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")
        handle.write(f"MOTIF {motif_name}\n")
        handle.write(
            f"letter-probability matrix: alength= 4 w= {len(pwm)} nsites= 1 E= 0\n"
        )
        for _, probs in pwm:
            handle.write(" ".join(f"{prob:.8f}" for prob in probs) + "\n")


def write_pfm_tsv(path: Path, pwm):
    with path.open("w") as handle:
        handle.write("position\tA\tC\tG\tT\n")
        for pos, probs in pwm:
            handle.write(
                f"{pos}\t{probs[0]:.10f}\t{probs[1]:.10f}\t{probs[2]:.10f}\t{probs[3]:.10f}\n"
            )


def check_pwm(pwm, start_pos: int, end_pos: int):
    expected_positions = list(range(start_pos, end_pos + 1))
    actual_positions = [pos for pos, _ in pwm]
    if actual_positions != expected_positions:
        raise ValueError(
            f"PWM positions {actual_positions} do not match expected {expected_positions}"
        )
    for pos, probs in pwm:
        total = sum(probs)
        if not math.isclose(total, 1.0, rel_tol=0.0, abs_tol=1e-8):
            raise ValueError(f"Position {pos} probabilities sum to {total}, expected 1.0")


def run_cmd(cmd, cwd=None):
    result = subprocess.run(
        cmd,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            + " ".join(str(part) for part in cmd)
            + "\nSTDOUT:\n"
            + result.stdout
            + "\nSTDERR:\n"
            + result.stderr
        )
    return result


def make_open_region_fasta(genome_fa: Path, open_region_bed: Path, fasta_out: Path):
    run_cmd(
        [
            str(BEDTOOLS),
            "getfasta",
            "-fi",
            str(genome_fa),
            "-bed",
            str(open_region_bed),
            "-fo",
            str(fasta_out),
        ]
    )


def run_fimo_scan(meme_file: Path, fasta_file: Path, fimo_dir: Path, thresh: float):
    if fimo_dir.exists():
        shutil.rmtree(fimo_dir)
    run_cmd(
        [
            str(FIMO),
            "--thresh",
            str(thresh),
            "--max-stored-scores",
            "100000000",
            "--verbosity",
            "5",
            "--oc",
            str(fimo_dir),
            str(meme_file),
            str(fasta_file),
        ]
    )


def parse_sequence_name(sequence_name: str):
    if ":" not in sequence_name:
        return sequence_name, None, None
    chrom, coords = sequence_name.split(":")
    start_s, end_s = coords.split("-")
    return chrom, int(start_s), int(end_s)


def fimo_to_bed(fimo_tsv: Path, motif_name: str, output_bed: Path):
    records = []
    with fimo_tsv.open() as handle:
        reader = csv.DictReader(
            (line for line in handle if line.strip() and not line.startswith("#")),
            delimiter="\t",
        )
        for row in reader:
            chrom, region_start, _ = parse_sequence_name(row["sequence_name"])
            if region_start is None:
                hit_start = int(row["start"]) - 1
                hit_end = int(row["stop"])
            else:
                hit_start = region_start + int(row["start"]) - 1
                hit_end = region_start + int(row["stop"])
            score = float(row["score"])
            records.append(
                (
                    chrom,
                    hit_start,
                    hit_end,
                    motif_name,
                    f"{score:.6f}",
                    row["strand"],
                )
            )

    if not records:
        return 0

    records.sort(key=lambda item: (item[0], item[1], item[2], item[5]))
    output_bed.parent.mkdir(parents=True, exist_ok=True)
    with output_bed.open("w") as handle:
        for record in records:
            handle.write("\t".join(str(field) for field in record) + "\n")
    return len(records)


def read_sample_sheet(sample_sheet: Path):
    sample_ids = []
    concentrations = []
    with sample_sheet.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample_ids.append(row["sample_id"])
            concentrations.append(row["concentration_M"])
    return sample_ids, concentrations


def write_top_tfbs_bed(input_bed: Path, output_bed: Path, top_n: int):
    records = []
    with input_bed.open() as handle:
        for order, line in enumerate(handle):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            try:
                score = float(fields[4])
            except ValueError:
                score = float("-inf")
            records.append((score, order, line))

    records.sort(key=lambda item: (-item[0], item[1]))
    chosen = records[:top_n]
    output_bed.parent.mkdir(parents=True, exist_ok=True)
    with output_bed.open("w") as handle:
        for _, _, line in chosen:
            handle.write(line)
    return len(records), len(chosen)


def run_plot(results_dir: Path, binding_sites: Path, sample_ids, concentrations, output_pdf: Path):
    sample_files = [
        str(results_dir / "01_compressed" / f"{sample_id}_sites.txt.gz")
        for sample_id in sample_ids
    ]
    missing = [path for path in sample_files if not Path(path).exists()]
    if missing:
        raise FileNotFoundError(f"Missing sample files: {missing}")

    pc_control = results_dir / "01_compressed" / "PC_control_sites.txt.gz"
    if not pc_control.exists():
        raise FileNotFoundError(f"Missing PC control file: {pc_control}")

    cmd = [
        str(PLOT_PYTHON),
        str(PLOT_SCRIPT),
        "--binding_sites",
        str(binding_sites),
        "--samples",
        ",".join(sample_files),
        "--sample_ids",
        ",".join(sample_ids),
        "--concentrations",
        ",".join(concentrations),
        "--pc_control",
        str(pc_control),
        "--output",
        str(output_pdf),
        "--motif_pick_ratio",
        "1.0",
        "--min_depth",
        "5",
        "--span",
        "100",
        "--plot_format",
        "pdf",
        "--plot_dpi",
        "300",
        "--normalize_flank",
    ]
    run_cmd(cmd)


def main():
    args = parse_args()

    results_tf = args.results_tf or args.tf_name
    results_dir = ANALYSIS_ROOT / "results" / results_tf
    config_path = args.config_path or (WORKFLOW_ROOT / "config" / f"config-{args.tf_name}.yaml")

    if args.binding_sites_bed is None and not args.logo_path.exists():
        raise FileNotFoundError(f"Missing logo TSV: {args.logo_path}")
    if not config_path.exists():
        raise FileNotFoundError(f"Missing config file: {config_path}")
    if not PLOT_PYTHON.exists():
        raise FileNotFoundError(f"Missing plotting python: {PLOT_PYTHON}")
    if not BEDTOOLS.exists():
        raise FileNotFoundError(f"Missing bedtools binary: {BEDTOOLS}")
    if not FIMO.exists():
        raise FileNotFoundError(f"Missing fimo binary: {FIMO}")

    config = read_config(config_path)
    genome_fa = Path(config["genomes"]["human"])
    open_region_bed = Path(config["open_regions"]["GM"])
    sample_sheet = WORKFLOW_ROOT / config["tfs"][args.tf_name]["sample_sheet"]

    if not genome_fa.exists():
        raise FileNotFoundError(f"Missing genome FASTA: {genome_fa}")
    if not open_region_bed.exists():
        raise FileNotFoundError(f"Missing GM open region BED: {open_region_bed}")
    if not sample_sheet.exists():
        raise FileNotFoundError(f"Missing sample sheet: {sample_sheet}")

    out_root = results_dir / args.output_tag
    motifs_dir = out_root / "05_motifs"
    binding_dir = out_root / "06_binding_sites"
    viz_dir = out_root / args.visualization_dir_name
    out_root.mkdir(parents=True, exist_ok=True)
    motifs_dir.mkdir(parents=True, exist_ok=True)
    binding_dir.mkdir(parents=True, exist_ok=True)
    viz_dir.mkdir(parents=True, exist_ok=True)

    motif_safe_tf = args.tf_name.replace("-", "_")
    motif_name = args.motif_name or (
        f"{motif_safe_tf}_linear_logo_pos{args.start_position}_{args.end_position}"
    )
    meme_file = motifs_dir / f"{motif_name}.meme"
    pfm_tsv = motifs_dir / f"{motif_name}.pfm.tsv"
    fasta_file = motifs_dir / "GM_open_regions.fa"
    fimo_dir = motifs_dir / f"{motif_name}.fimo_out"
    fimo_bed = motifs_dir / f"{motif_name}.bed"
    binding_sites = binding_dir / f"{results_tf}.{args.output_tag}.open_region.bed"
    output_pdf_suffix = (
        f".top{args.top_n_tfbs}TFBS" if args.top_n_tfbs is not None else ""
    )
    output_pdf = viz_dir / f"{results_tf}.{args.output_tag}{output_pdf_suffix}.conversion_plot.pdf"
    plot_binding_sites = binding_sites
    top_tfbs_bed = (
        viz_dir / f"{results_tf}.{args.output_tag}.top{args.top_n_tfbs}TFBS.open_region.bed"
        if args.top_n_tfbs is not None
        else None
    )
    summary_tsv = out_root / f"summary{output_pdf_suffix}.tsv"

    if args.binding_sites_bed is None:
        logo_rows = read_linear_logo(args.logo_path, args.start_position, args.end_position)
        pwm = pwm_from_logo(logo_rows)
        check_pwm(pwm, args.start_position, args.end_position)
        write_meme(meme_file, motif_name, pwm)
        write_pfm_tsv(pfm_tsv, pwm)

        make_open_region_fasta(genome_fa, open_region_bed, fasta_file)
        run_fimo_scan(meme_file, fasta_file, fimo_dir, args.fimo_thresh)

        fimo_tsv = fimo_dir / "fimo.tsv"
        if not fimo_tsv.exists():
            raise FileNotFoundError(f"FIMO did not produce output TSV: {fimo_tsv}")

        hit_count = fimo_to_bed(fimo_tsv, motif_name, fimo_bed)
        if hit_count == 0:
            raise RuntimeError(
                f"FIMO returned 0 hits at threshold {args.fimo_thresh}. No conversion plot was generated."
            )

        shutil.copyfile(fimo_bed, binding_sites)
    else:
        if not args.binding_sites_bed.exists():
            raise FileNotFoundError(f"Missing binding sites BED: {args.binding_sites_bed}")
        if args.binding_sites_bed.resolve() != binding_sites.resolve():
            shutil.copyfile(args.binding_sites_bed, binding_sites)
        hit_count = sum(1 for line in binding_sites.open() if line.strip() and not line.startswith("#"))

    plot_hit_count = hit_count
    if args.top_n_tfbs is not None:
        if args.top_n_tfbs <= 0:
            raise ValueError("--top-n-tfbs must be positive")
        hit_count, plot_hit_count = write_top_tfbs_bed(
            binding_sites, top_tfbs_bed, args.top_n_tfbs
        )
        if plot_hit_count == 0:
            raise RuntimeError(f"No TFBS records available in {binding_sites}")
        plot_binding_sites = top_tfbs_bed

    sample_ids, concentrations = read_sample_sheet(sample_sheet)
    run_plot(results_dir, plot_binding_sites, sample_ids, concentrations, output_pdf)

    if not output_pdf.exists() or output_pdf.stat().st_size == 0:
        raise RuntimeError(f"Plot generation failed: {output_pdf}")

    with summary_tsv.open("w") as handle:
        handle.write("key\tvalue\n")
        handle.write(f"tf_name\t{args.tf_name}\n")
        handle.write(f"results_tf\t{results_tf}\n")
        handle.write(f"config_path\t{config_path}\n")
        handle.write(f"logo_tsv\t{args.logo_path if args.binding_sites_bed is None else ''}\n")
        handle.write(f"positions\t{args.start_position}-{args.end_position}\n")
        handle.write(f"meme_file\t{meme_file if args.binding_sites_bed is None else ''}\n")
        handle.write(f"pfm_tsv\t{pfm_tsv if args.binding_sites_bed is None else ''}\n")
        handle.write(f"open_region_bed\t{open_region_bed}\n")
        handle.write(f"open_region_fasta\t{fasta_file if args.binding_sites_bed is None else ''}\n")
        handle.write(f"fimo_dir\t{fimo_dir if args.binding_sites_bed is None else ''}\n")
        handle.write(f"fimo_bed\t{fimo_bed if args.binding_sites_bed is None else ''}\n")
        handle.write(f"source_binding_sites_bed\t{args.binding_sites_bed or ''}\n")
        handle.write(f"binding_sites_bed\t{binding_sites}\n")
        handle.write(f"plot_binding_sites_bed\t{plot_binding_sites}\n")
        handle.write(f"conversion_plot_pdf\t{output_pdf}\n")
        handle.write(f"fimo_thresh\t{args.fimo_thresh}\n")
        handle.write(f"hit_count\t{hit_count}\n")
        handle.write(f"top_n_tfbs\t{args.top_n_tfbs or ''}\n")
        handle.write(f"plot_hit_count\t{plot_hit_count}\n")

    print(f"[OK] MEME motif: {meme_file}")
    print(f"[OK] FIMO hits: {hit_count}")
    print(f"[OK] Binding sites BED: {binding_sites}")
    print(f"[OK] Plot binding sites BED: {plot_binding_sites}")
    print(f"[OK] Conversion plot: {output_pdf}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

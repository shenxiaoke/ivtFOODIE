#!/usr/bin/env python3

from __future__ import annotations
import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


def find_motif_bed(motif_dir: Path, tf_name: str) -> Path | None:
    """
    Find corresponding motif BED file.

    Examples:
        MLX     -> MLX.bed
        MLX-FL  -> MLX.bed
        SP1-FL  -> SP1.bed

    Matching is case-insensitive.
    """
    # base_name = tf_name.replace("-FL", "")
    base_name = tf_name.split('-')[0]

    motif_files = list(motif_dir.glob("*.bed"))
    name_to_path = {p.stem.lower(): p for p in motif_files}

    return name_to_path.get(base_name.lower(), None)


def read_motif_bed_with_uid(motif_bed: Path) -> tuple[pd.DataFrame, list[str]]:
    """
    Read motif BED and add a unique motif_uid column.
    """
    motif_df = pd.read_csv(
        motif_bed,
        sep="\t",
        header=None,
        dtype={0: str},
        comment="#"
    )

    n_cols = motif_df.shape[1]
    motif_cols = [f"motif_col{i+1}" for i in range(n_cols)]

    motif_df.columns = motif_cols
    motif_df["motif_uid"] = np.arange(len(motif_df), dtype=np.int64)

    return motif_df, motif_cols


def prepare_site_bed(fitting_tsv: Path) -> pd.DataFrame:
    """
    Read single C/G fitting Kd table.

    Required columns:
        ch, st, ed, Kd, spearmanr

    Output BED-like dataframe:
        ch, st, ed, site_uid, neglog10Kd, spearmanr, Kd
    """
    df = pd.read_csv(fitting_tsv, sep="\t")

    required = ["ch", "st", "ed", "Kd", "spearmanr"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{fitting_tsv} is missing required columns: {missing}"
        )

    df = df[required].copy()

    df["Kd"] = pd.to_numeric(df["Kd"], errors="coerce")
    df["spearmanr"] = pd.to_numeric(df["spearmanr"], errors="coerce")
    df["st"] = pd.to_numeric(df["st"], errors="coerce")
    df["ed"] = pd.to_numeric(df["ed"], errors="coerce")

    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["ch", "st", "ed", "Kd"])

    # Kd must be positive before log transform
    df = df[df["Kd"] > 0].copy()

    df["st"] = df["st"].astype(np.int64)
    df["ed"] = df["ed"].astype(np.int64)

    df["site_uid"] = np.arange(len(df), dtype=np.int64)
    df["neglog10Kd"] = -np.log10(df["Kd"])

    site_bed = df[
        ["ch", "st", "ed", "site_uid", "neglog10Kd", "spearmanr", "Kd"]
    ].copy()

    return site_bed


def run_bedtools_intersect(
    motif_df: pd.DataFrame,
    motif_cols: list[str],
    site_bed: pd.DataFrame,
    tmpdir: Path,
) -> pd.DataFrame:
    """
    Run:
        bedtools intersect -a motif_with_uid.bed -b site.bed -wa -wb
    """
    motif_tmp = tmpdir / "motif.with_uid.bed"
    site_tmp = tmpdir / "sites.kd.bed"
    out_tmp = tmpdir / "intersect.tsv"

    motif_df.to_csv(
        motif_tmp,
        sep="\t",
        header=False,
        index=False,
    )

    site_bed.to_csv(
        site_tmp,
        sep="\t",
        header=False,
        index=False,
    )

    cmd = [
        "bedtools",
        "intersect",
        "-a", str(motif_tmp),
        "-b", str(site_tmp),
        "-wa",
        "-wb",
    ]

    with open(out_tmp, "w") as fout:
        subprocess.run(cmd, check=True, stdout=fout)

    a_cols = motif_cols + ["motif_uid"]
    b_cols = [
        "site_ch",
        "site_st",
        "site_ed",
        "site_uid",
        "site_neglog10Kd",
        "site_spearmanr",
        "site_Kd",
    ]

    if out_tmp.stat().st_size == 0:
        return pd.DataFrame(columns=a_cols + b_cols)

    inter = pd.read_csv(
        out_tmp,
        sep="\t",
        header=None,
        names=a_cols + b_cols,
    )

    return inter


def summarize_to_motif_level(
    motif_df: pd.DataFrame,
    motif_cols: list[str],
    inter: pd.DataFrame,
    keep_all: bool = False,
) -> pd.DataFrame:
    """
    Summarize C/G-level Kd to motif-level Kd.
    """
    motif_meta = motif_df[motif_cols + ["motif_uid"]].copy()

    if inter.empty:
        out = motif_meta.copy()
        out["n_CG"] = 0
        out["motif_neglog10Kd_mean"] = np.nan
        out["motif_Kd_from_mean_neglog10Kd"] = np.nan
        out["motif_spearmanr_mean"] = np.nan
        return out if keep_all else out.iloc[0:0].copy()

    inter["site_neglog10Kd"] = pd.to_numeric(
        inter["site_neglog10Kd"], errors="coerce"
    )
    inter["site_spearmanr"] = pd.to_numeric(
        inter["site_spearmanr"], errors="coerce"
    )
    inter["site_Kd"] = pd.to_numeric(
        inter["site_Kd"], errors="coerce"
    )

    agg = (
        inter
        .groupby("motif_uid", as_index=False)
        .agg(
            n_CG=("site_uid", "count"),
            motif_neglog10Kd_mean=("site_neglog10Kd", "mean"),
            motif_spearmanr_mean=("site_spearmanr", "mean"),
            site_Kd_direct_mean=("site_Kd", "mean"),
            site_neglog10Kd_min=("site_neglog10Kd", "min"),
            site_neglog10Kd_max=("site_neglog10Kd", "max"),
        )
    )

    agg["motif_Kd_from_mean_neglog10Kd"] = (
        10 ** (-agg["motif_neglog10Kd_mean"])
    )

    if keep_all:
        out = motif_meta.merge(agg, on="motif_uid", how="left")
        out["n_CG"] = out["n_CG"].fillna(0).astype(int)
    else:
        out = motif_meta.merge(agg, on="motif_uid", how="inner")

    # Put important summary columns earlier
    summary_cols = [
        "n_CG",
        "motif_neglog10Kd_mean",
        "motif_Kd_from_mean_neglog10Kd",
        "motif_spearmanr_mean",
        "site_Kd_direct_mean",
        "site_neglog10Kd_min",
        "site_neglog10Kd_max",
    ]

    out = out[motif_cols + ["motif_uid"] + summary_cols]

    return out


def process_one_fitting_file(
    fitting_tsv: Path,
    motif_bed: Path,
    output_tsv: Path,
    keep_all: bool = False,
):
    motif_df, motif_cols = read_motif_bed_with_uid(motif_bed)
    site_bed = prepare_site_bed(fitting_tsv)

    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)
        inter = run_bedtools_intersect(
            motif_df=motif_df,
            motif_cols=motif_cols,
            site_bed=site_bed,
            tmpdir=tmpdir,
        )

    motif_level = summarize_to_motif_level(
        motif_df=motif_df,
        motif_cols=motif_cols,
        inter=inter,
        keep_all=keep_all,
    )

    motif_level.to_csv(output_tsv, sep="\t", index=False)

    print(
        f"[OK] {fitting_tsv} -> {output_tsv} "
        f"| motif rows: {len(motif_level)}"
    )


def process_one_tf(
    tf_dir: Path,
    motif_dir: Path,
    motif_tf_name: str | None = None,
    keep_all: bool = False,
):
    tf_name = tf_dir.name
    motif_lookup_name = motif_tf_name or tf_name
    motif_bed = find_motif_bed(motif_dir, motif_lookup_name)

    if motif_bed is None:
        print(
            f"[SKIP] No motif BED found for {tf_name} "
            f"(motif TF name: {motif_lookup_name})"
        )
        return

    fitting_dir = tf_dir / "03_fitting"
    out_dir = tf_dir / "08_fitting_merge"

    if not fitting_dir.exists():
        print(f"[SKIP] No 03_fitting dir: {tf_dir}")
        return

    fitting_files = [
        (
            fitting_dir / f"{tf_name}.fittingKd.basic.tsv",
            out_dir / f"{tf_name}.motif_level.fittingKd.basic.tsv",
        ),
        (
            fitting_dir / f"{tf_name}.fittingKd.flank_corrected.tsv",
            out_dir / f"{tf_name}.motif_level.fittingKd.flank_corrected.tsv",
        ),
    ]

    print(f"\n[TF] {tf_name}")
    print(f"     motif TF name: {motif_lookup_name}")
    print(f"     motif: {motif_bed}")

    for fitting_tsv, output_tsv in fitting_files:
        if not fitting_tsv.exists():
            print(f"[SKIP] Missing fitting file: {fitting_tsv}")
            continue

        process_one_fitting_file(
            fitting_tsv=fitting_tsv,
            motif_bed=motif_bed,
            output_tsv=output_tsv,
            keep_all=keep_all,
        )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Merge single C/G-level Kd into motif-level Kd. "
            "Kd is first transformed to -log10(Kd), then averaged within motif."
        )
    )

    parser.add_argument(
        "--results-dir",
        required=True,
        help="results directory containing TF subdirectories, e.g. .../results",
    )

    parser.add_argument(
        "--motif-dir",
        required=True,
        help="directory containing split JASPAR BED files, e.g. JASPAR_split_d50",
    )

    parser.add_argument(
        "--tf",
        nargs="*",
        default=None,
        help=(
            "Case-insensitive keywords used to match TF directory names. "
            "If not provided, process all subdirectories under results-dir."
        ),
    )

    parser.add_argument(
        "--tf-name",
        "--tf_name",
        "--TF_name",
        dest="tf_name",
        default=None,
        help=(
            "TF name used to select the motif BED (for example SP1, YY1, or "
            "CTCF). When provided, the same motif BED is used for every "
            "directory matched by --tf."
        ),
    )

    parser.add_argument(
        "--keep-all",
        action="store_true",
        help=(
            "Keep motif rows even if no C/G fitting site overlaps. "
            "Default: only output motifs with at least one overlapping C/G."
        ),
    )

    args = parser.parse_args()

    results_dir = Path(args.results_dir).resolve()
    motif_dir = Path(args.motif_dir).resolve()

    if args.tf_name and not args.tf:
        parser.error("--tf-name requires at least one --tf keyword")

    candidate_tf_dirs = [
        p for p in results_dir.iterdir()
        if p.is_dir() and (p / "03_fitting").exists()
    ]

    if args.tf:
        tf_keywords = [keyword.lower() for keyword in args.tf]
        tf_dirs = [
            p for p in candidate_tf_dirs
            if any(keyword in p.name.lower() for keyword in tf_keywords)
        ]

        if not tf_dirs:
            parser.error(
                "No TF directories containing 03_fitting matched --tf "
                f"keywords: {args.tf}"
            )
    else:
        tf_dirs = candidate_tf_dirs

    for tf_dir in sorted(tf_dirs):
        process_one_tf(
            tf_dir=tf_dir,
            motif_dir=motif_dir,
            motif_tf_name=args.tf_name,
            keep_all=args.keep_all,
        )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot TF motif-centered conversion profiles.

For every TF-matched sample under ratio_dir, this script reads:
    <ratio_dir>/<sample>/<sample>_sites.txt.gz

It computes two kinds of profiles:
    1. Positive-control-corrected profile:
         sample conversion ratio / positive-control conversion ratio
    2. Uncorrected profile:
         sample conversion ratio only

The output PDF contains two pages:
    Page 1: PC-corrected profiles
        - top subplot: no further normalization
        - bottom subplot: flank-normalized
    Page 2: uncorrected sample profiles
        - top subplot: no further normalization
        - bottom subplot: flank-normalized

A TSV containing all plotted values is also written for reproducibility.
"""

from __future__ import annotations

import argparse
import gzip
import math
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

try:
    import tabix
except ImportError as exc:
    raise SystemExit(
        "Cannot import 'tabix'. Install pytabix in the current environment, for example:\n"
        "    conda install -c bioconda pytabix\n"
        "or:\n"
        "    pip install pytabix"
    ) from exc


SCRIPT_VERSION = "2026-07-13-v4-output-tag"


DEFAULT_MOTIF_DIR = Path(
    "/dshare/xielab/analysisdata/TFC/shenk/important_data/TFmotifs/"
    "jaspar/JASPAR_split_d50"
)


def eprint(*args, **kwargs) -> None:
    print(*args, file=sys.stderr, **kwargs)


def natural_key(text: str) -> List[object]:
    """Natural-sort key: 200nM appears before 1000nM."""
    return [
        int(piece) if piece.isdigit() else piece.lower()
        for piece in re.split(r"(\d+)", text)
    ]


def sanitize_filename(text: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", text.strip())
    return cleaned.strip("._-") or "output"


def resolve_sample_filter_label(
    sample_regex: Optional[str],
    sample_regex_name: Optional[str],
) -> Optional[str]:
    """
    Return the label used to distinguish SAMPLE_REGEX-specific outputs.

    Priority:
      1. Explicit --sample-regex-name / SAMPLE_REGEX_NAME.
      2. A filename-safe form of --sample-regex / SAMPLE_REGEX.
      3. None when no custom regex is used.
    """
    if sample_regex_name is not None:
        label = sample_regex_name.strip()
        if not label:
            raise ValueError("sample_regex_name cannot be empty")
        return label
    if sample_regex:
        return sanitize_filename(sample_regex)
    return None


def open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def resolve_motif_bed(tf_name: str, motif_bed: Optional[str], motif_dir: Path) -> Path:
    """
    Resolve the motif BED file.

    If --motif-bed is absent, match the file stem to TF name case-insensitively.
    Therefore NRF1 can resolve to Nrf1.bed, while CTCF will not accidentally
    resolve to CTCFL.bed.
    """
    if motif_bed:
        path = Path(motif_bed).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Motif BED does not exist: {path}")
        return path

    motif_dir = motif_dir.expanduser().resolve()
    if not motif_dir.is_dir():
        raise FileNotFoundError(f"Motif directory does not exist: {motif_dir}")

    tf_folded = tf_name.casefold()
    candidates: List[Path] = []
    for path in motif_dir.iterdir():
        if not path.is_file():
            continue
        lower_name = path.name.casefold()
        if lower_name.endswith(".bed"):
            stem = path.name[:-4]
        elif lower_name.endswith(".bed.gz"):
            stem = path.name[:-7]
        else:
            continue
        if stem.casefold() == tf_folded:
            candidates.append(path.resolve())

    candidates.sort(key=lambda p: natural_key(p.name))
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise FileNotFoundError(
            f"No case-insensitive exact motif BED match for TF '{tf_name}' in {motif_dir}.\n"
            "Provide the BED explicitly as the third argument of the bash script."
        )
    raise RuntimeError(
        "Multiple motif BED files matched the TF name exactly:\n  "
        + "\n  ".join(str(p) for p in candidates)
    )


def resolve_sample_dir(ratio_dir: Path, sample_name: str) -> Path:
    """Resolve a sample directory by exact, case-insensitive, then unique substring match."""
    exact = ratio_dir / sample_name
    if exact.is_dir():
        return exact.resolve()

    dirs = [p for p in ratio_dir.iterdir() if p.is_dir()]
    ci_exact = [p for p in dirs if p.name.casefold() == sample_name.casefold()]
    if len(ci_exact) == 1:
        return ci_exact[0].resolve()

    substring = [p for p in dirs if sample_name.casefold() in p.name.casefold()]
    if len(substring) == 1:
        return substring[0].resolve()
    if not substring:
        raise FileNotFoundError(
            f"Positive-control sample directory was not found under {ratio_dir}: {sample_name}"
        )
    raise RuntimeError(
        f"Positive-control name '{sample_name}' matched multiple directories. "
        "Please provide the complete directory name:\n  "
        + "\n  ".join(sorted((p.name for p in substring), key=natural_key))
    )


def resolve_sites_file(sample_dir: Path) -> Path:
    """Resolve the unique bgzip/tabix-indexed *_sites.txt.gz file for one sample."""
    preferred = sample_dir / f"{sample_dir.name}_sites.txt.gz"
    if preferred.is_file():
        sites = preferred
    else:
        candidates = sorted(sample_dir.glob("*_sites.txt.gz"), key=lambda p: natural_key(p.name))
        if len(candidates) != 1:
            raise RuntimeError(
                f"Expected exactly one *_sites.txt.gz under {sample_dir}, found {len(candidates)}:\n  "
                + "\n  ".join(str(p) for p in candidates)
            )
        sites = candidates[0]

    tbi = Path(str(sites) + ".tbi")
    if not tbi.is_file():
        raise FileNotFoundError(f"Tabix index is missing: {tbi}")
    return sites.resolve()


def build_tf_sample_pattern(tf_name: str) -> re.Pattern[str]:
    """
    Match TF as a sample-name token, optionally followed by FL or DBD.

    Examples for TF=SP1:
      matches: ...-SP1-..., ...-SP1FL-..., ...-SP1DBD-...
      does not match: ...-SP10-...

    Examples for TF=CTCF:
      matches CTCF and CTCFDBD, but does not match CTCFL.
    """
    escaped = re.escape(tf_name)
    return re.compile(
        rf"(?<![A-Za-z0-9]){escaped}(?:FL|DBD)?(?=[^A-Za-z0-9]|$)",
        flags=re.IGNORECASE,
    )


def discover_tf_samples(
    ratio_dir: Path,
    tf_name: str,
    control_dir: Path,
    sample_regex: Optional[str] = None,
) -> List[Tuple[str, Path]]:
    pattern = re.compile(sample_regex, re.IGNORECASE) if sample_regex else build_tf_sample_pattern(tf_name)
    matched: List[Tuple[str, Path]] = []

    for sample_dir in ratio_dir.iterdir():
        if not sample_dir.is_dir() or sample_dir.resolve() == control_dir.resolve():
            continue
        if pattern.search(sample_dir.name) is None:
            continue
        try:
            sites = resolve_sites_file(sample_dir)
        except Exception as exc:
            eprint(f"[WARN] Skipping {sample_dir.name}: {exc}")
            continue
        matched.append((sample_dir.name, sites))

    matched.sort(key=lambda item: natural_key(item[0]))
    if not matched:
        regex_note = f" using regex {sample_regex!r}" if sample_regex else ""
        raise FileNotFoundError(
            f"No usable TF sample directories matching '{tf_name}' were found under "
            f"{ratio_dir}{regex_note}."
        )
    return matched


def load_and_select_motifs(
    motif_bed: Path,
    motif_pick_ratio: float,
    pick: str,
) -> Tuple[List[Tuple[str, int, int, str, float, str]], float, int]:
    if not (0 < motif_pick_ratio <= 1):
        raise ValueError("motif_pick_ratio must satisfy 0 < ratio <= 1")

    motifs: List[Tuple[str, int, int, str, float, str]] = []
    scores: List[float] = []
    malformed = 0

    with open_text(motif_bed) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                malformed += 1
                continue
            try:
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                name = fields[3]
                score = float(fields[4])
                strand = fields[5]
            except ValueError:
                malformed += 1
                continue
            if end <= start or strand not in {"+", "-"} or not math.isfinite(score):
                malformed += 1
                continue
            motifs.append((chrom, start, end, name, score, strand))
            scores.append(score)

    if not motifs:
        raise RuntimeError(f"No valid motif records were read from {motif_bed}")

    score_array = np.asarray(scores, dtype=float)
    if pick == "top":
        threshold = float(np.percentile(score_array, 100.0 * (1.0 - motif_pick_ratio)))
        selected = [m for m in motifs if m[4] >= threshold]
    else:
        threshold = float(np.percentile(score_array, 100.0 * motif_pick_ratio))
        selected = [m for m in motifs if m[4] <= threshold]

    return selected, threshold, malformed


def parse_site_records(records: Iterable[Sequence[str]]) -> Dict[int, Tuple[int, float]]:
    """
    Convert tabix records to {position: (depth, conversion_ratio)}.

    Expected columns follow the user's existing pipeline:
      record[1] = genomic position
      record[3] + record[4] = depth
      record[5] = conversion ratio
    """
    parsed: Dict[int, Tuple[int, float]] = {}
    for record in records:
        if len(record) < 6:
            continue
        try:
            pos = int(record[1])
            depth = int(record[3]) + int(record[4])
            ratio = float(record[5])
        except (TypeError, ValueError):
            continue
        if math.isfinite(ratio):
            parsed[pos] = (depth, ratio)
    return parsed


def query_sites(handle, chrom: str, start: int, end: int) -> Dict[int, Tuple[int, float]]:
    start = max(0, int(start))
    end = max(start + 1, int(end))
    try:
        return parse_site_records(handle.query(chrom, start, end))
    except Exception:
        return {}


def aggregate(values: List[float], method: str) -> float:
    if not values:
        return float("nan")
    array = np.asarray(values, dtype=float)
    if method == "median":
        return float(np.nanmedian(array))
    return float(np.nanmean(array))


def compute_profiles(
    sample_files: Sequence[Tuple[str, Path]],
    control_file: Path,
    motifs: Sequence[Tuple[str, int, int, str, float, str]],
    span: int,
    min_depth: int,
    method: str,
) -> Tuple[np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, int], Dict[str, int]]:
    """
    Compute both PC-corrected and uncorrected profiles.

    Returns:
        x
        corrected_profiles
        raw_profiles
        corrected_contribution_counts
        raw_contribution_counts
    """
    if span <= 0 or span % 2 != 0:
        raise ValueError("span must be a positive even integer")
    if min_depth < 0:
        raise ValueError("min_depth must be >= 0")

    half = span // 2
    x = np.arange(-half, half, dtype=int)

    corrected_values: Dict[str, List[List[float]]] = {
        sample_name: [[] for _ in range(span)] for sample_name, _ in sample_files
    }
    raw_values: Dict[str, List[List[float]]] = {
        sample_name: [[] for _ in range(span)] for sample_name, _ in sample_files
    }
    corrected_contribution_counts = {sample_name: 0 for sample_name, _ in sample_files}
    raw_contribution_counts = {sample_name: 0 for sample_name, _ in sample_files}

    control_handle = tabix.open(str(control_file))
    sample_handles = {
        sample_name: tabix.open(str(sites_file))
        for sample_name, sites_file in sample_files
    }

    total = len(motifs)
    for motif_index, (chrom, start, end, _name, _score, strand) in enumerate(motifs, start=1):
        if motif_index == 1 or motif_index % 1000 == 0 or motif_index == total:
            print(f"[INFO] Processing motif {motif_index:,}/{total:,}", end="\r", flush=True)

        center = (start + end) // 2
        query_start = center - half
        query_end = center + half
        control_records = query_sites(control_handle, chrom, query_start, query_end)

        for sample_name, sample_handle in sample_handles.items():
            sample_records = query_sites(sample_handle, chrom, query_start, query_end)
            if not sample_records:
                continue

            for pos, (sample_depth, sample_ratio) in sample_records.items():
                if sample_depth < min_depth:
                    continue
                if not math.isfinite(sample_ratio):
                    continue

                if strand == "+":
                    index = pos - center + half
                else:
                    index = center - pos + half

                if not (0 <= index < span):
                    continue

                raw_values[sample_name][index].append(sample_ratio)
                raw_contribution_counts[sample_name] += 1

                if pos not in control_records:
                    continue
                control_depth, control_ratio = control_records[pos]
                if control_depth < min_depth:
                    continue
                if control_ratio == 0 or not math.isfinite(control_ratio):
                    continue

                corrected = sample_ratio / control_ratio
                if not math.isfinite(corrected):
                    continue

                corrected_values[sample_name][index].append(corrected)
                corrected_contribution_counts[sample_name] += 1

    print(" " * 100, end="\r")

    corrected_profiles = {
        sample_name: np.asarray(
            [aggregate(position_values, method) for position_values in per_position],
            dtype=float,
        )
        for sample_name, per_position in corrected_values.items()
    }
    raw_profiles = {
        sample_name: np.asarray(
            [aggregate(position_values, method) for position_values in per_position],
            dtype=float,
        )
        for sample_name, per_position in raw_values.items()
    }
    return x, corrected_profiles, raw_profiles, corrected_contribution_counts, raw_contribution_counts


def flank_normalize(
    x: np.ndarray,
    profile: np.ndarray,
    flank_exclusion: int,
) -> Tuple[np.ndarray, float]:
    """Normalize by the mean outside [-flank_exclusion, flank_exclusion)."""
    flank_mask = (x < -flank_exclusion) | (x >= flank_exclusion)
    finite_mask = flank_mask & np.isfinite(profile)
    if not np.any(finite_mask):
        return np.full_like(profile, np.nan, dtype=float), float("nan")

    factor = float(np.nanmean(profile[finite_mask]))
    if not math.isfinite(factor) or factor == 0:
        return np.full_like(profile, np.nan, dtype=float), factor
    return profile / factor, factor


def compact_plot_labels(sample_names: Sequence[str]) -> Dict[str, str]:
    labels: Dict[str, str] = {}
    used: Dict[str, int] = {}

    for sample_name in sample_names:
        label = re.sub(r"^.*?ivtfoodie[-_]", "", sample_name, flags=re.IGNORECASE)
        label = re.sub(r"_L\d+_[A-Za-z0-9]+$", "", label)
        label = re.sub(r"_L\d+$", "", label)
        label = label or sample_name

        used[label] = used.get(label, 0) + 1
        if used[label] > 1:
            label = f"{label} [{used[label]}]"
        labels[sample_name] = label
    return labels


def write_profile_table(
    output_tsv: Path,
    x: np.ndarray,
    corrected_profiles: Dict[str, np.ndarray],
    corrected_normalized_profiles: Dict[str, np.ndarray],
    raw_profiles: Dict[str, np.ndarray],
    raw_normalized_profiles: Dict[str, np.ndarray],
) -> None:
    sample_names = list(corrected_profiles)
    with output_tsv.open("w") as handle:
        columns = ["relative_position"]
        for name in sample_names:
            columns.extend(
                [
                    f"{name}__pc_corrected",
                    f"{name}__pc_corrected_flank_normalized",
                    f"{name}__raw",
                    f"{name}__raw_flank_normalized",
                ]
            )
        handle.write("\t".join(columns) + "\n")

        for i, pos in enumerate(x):
            row = [str(int(pos))]
            for name in sample_names:
                values = [
                    corrected_profiles[name][i],
                    corrected_normalized_profiles[name][i],
                    raw_profiles[name][i],
                    raw_normalized_profiles[name][i],
                ]
                for value in values:
                    row.append("nan" if not np.isfinite(value) else f"{value:.10g}")
            handle.write("\t".join(row) + "\n")


def _draw_two_subplots_page(
    pdf: PdfPages,
    tf_name: str,
    motif_bed: Path,
    control_name: str,
    sample_filter_label: Optional[str],
    sample_regex: Optional[str],
    x: np.ndarray,
    profiles: Dict[str, np.ndarray],
    normalized_profiles: Dict[str, np.ndarray],
    sample_names: Sequence[str],
    page_title: str,
    raw_title: str,
    norm_title: str,
    y_label: str,
    motif_count: int,
    score_threshold: float,
    motif_pick_ratio: float,
    pick: str,
    flank_exclusion: int,
    raw_ylim: Optional[Tuple[float, float]],
    normalized_ylim: Optional[Tuple[float, float]],
) -> None:
    labels = compact_plot_labels(sample_names)

    # Previous height was 10. Increase to 1.5x -> 15.
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(13, 15), sharex=True)
    ax_raw, ax_norm = axes

    for sample_name in sample_names:
        ax_raw.plot(x, profiles[sample_name], linewidth=1.6, label=labels[sample_name])
        ax_norm.plot(
            x,
            normalized_profiles[sample_name],
            linewidth=1.6,
            label=labels[sample_name],
        )

    for ax in axes:
        ax.axvline(0, linestyle="--", linewidth=1.0, alpha=0.7)
        ax.axhline(1, linestyle=":", linewidth=1.0, alpha=0.7)
        ax.grid(alpha=0.2)
        ax.set_ylabel(y_label)

    ax_raw.set_title(raw_title)
    ax_norm.set_title(norm_title)
    ax_norm.set_xlabel("Position relative to motif center (bp)")

    if raw_ylim is not None:
        ax_raw.set_ylim(raw_ylim)
    if normalized_ylim is not None:
        ax_norm.set_ylim(normalized_ylim)

    handles, legend_labels = ax_raw.get_legend_handles_labels()
    fig.legend(
        handles,
        legend_labels,
        loc="center left",
        bbox_to_anchor=(0.82, 0.5),
        frameon=False,
        fontsize=8,
    )

    selection_text = (
        f"{pick} {motif_pick_ratio:.1%}; score threshold={score_threshold:g}; "
        f"selected motifs={motif_count:,}"
    )
    filter_text = ""
    if sample_filter_label is not None:
        filter_text = f"\nSample filter: {sample_filter_label}"
        if sample_regex and sample_regex != sample_filter_label:
            filter_text += f" (regex: {sample_regex})"

    fig.suptitle(
        f"{tf_name} motif-centered conversion\n"
        f"{page_title}\n"
        f"Control: {control_name}{filter_text}\n"
        f"Motif BED: {motif_bed.name}; {selection_text}",
        fontsize=13,
    )
    fig.subplots_adjust(left=0.09, right=0.80, top=0.90, bottom=0.07, hspace=0.32)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def plot_profiles(
    output_pdf: Path,
    tf_name: str,
    motif_bed: Path,
    control_name: str,
    sample_filter_label: Optional[str],
    sample_regex: Optional[str],
    x: np.ndarray,
    corrected_profiles: Dict[str, np.ndarray],
    corrected_normalized_profiles: Dict[str, np.ndarray],
    raw_profiles: Dict[str, np.ndarray],
    raw_normalized_profiles: Dict[str, np.ndarray],
    motif_count: int,
    score_threshold: float,
    motif_pick_ratio: float,
    pick: str,
    flank_exclusion: int,
    raw_ylim: Optional[Tuple[float, float]],
    normalized_ylim: Optional[Tuple[float, float]],
) -> None:
    sample_names = list(corrected_profiles)
    with PdfPages(output_pdf) as pdf:
        _draw_two_subplots_page(
            pdf=pdf,
            tf_name=tf_name,
            motif_bed=motif_bed,
            control_name=control_name,
            sample_filter_label=sample_filter_label,
            sample_regex=sample_regex,
            x=x,
            profiles=corrected_profiles,
            normalized_profiles=corrected_normalized_profiles,
            sample_names=sample_names,
            page_title="Page 1: positive-control-corrected profiles",
            raw_title="PC-corrected conversion profile (no further normalization)",
            norm_title=(
                f"PC-corrected profile normalized by flank mean "
                f"(outside center ±{flank_exclusion} bp)"
            ),
            y_label="Conversion ratio",
            motif_count=motif_count,
            score_threshold=score_threshold,
            motif_pick_ratio=motif_pick_ratio,
            pick=pick,
            flank_exclusion=flank_exclusion,
            raw_ylim=raw_ylim,
            normalized_ylim=normalized_ylim,
        )

        _draw_two_subplots_page(
            pdf=pdf,
            tf_name=tf_name,
            motif_bed=motif_bed,
            control_name=control_name,
            sample_filter_label=sample_filter_label,
            sample_regex=sample_regex,
            x=x,
            profiles=raw_profiles,
            normalized_profiles=raw_normalized_profiles,
            sample_names=sample_names,
            page_title="Page 2: uncorrected sample profiles",
            raw_title="Raw conversion profile (no positive-control correction)",
            norm_title=(
                f"Raw profile normalized by flank mean "
                f"(outside center ±{flank_exclusion} bp)"
            ),
            y_label="Conversion ratio",
            motif_count=motif_count,
            score_threshold=score_threshold,
            motif_pick_ratio=motif_pick_ratio,
            pick=pick,
            flank_exclusion=flank_exclusion,
            raw_ylim=raw_ylim,
            normalized_ylim=normalized_ylim,
        )


def parse_ylim(values: Optional[Sequence[float]]) -> Optional[Tuple[float, float]]:
    if values is None:
        return None
    low, high = float(values[0]), float(values[1])
    if not low < high:
        raise ValueError("ylim requires LOW < HIGH")
    return low, high


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Find all ratio-folder samples matching a TF, calculate motif-centered "
            "conversion profiles, and write a two-page PDF. "
            "Page 1 shows PC-corrected profiles; page 2 shows raw uncorrected profiles."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--tf", required=True, help="TF name, for example CTCF, NRF1, SP1")
    parser.add_argument(
        "--positive-control",
        required=True,
        help="Positive-control sample directory name under ratio_dir",
    )
    parser.add_argument(
        "--motif-bed",
        default=None,
        help="Explicit TFBS BED or BED.GZ path. If omitted, resolve from --motif-dir.",
    )
    parser.add_argument(
        "--motif-dir",
        type=Path,
        default=DEFAULT_MOTIF_DIR,
        help="Default motif BED directory",
    )
    parser.add_argument(
        "--ratio-dir",
        type=Path,
        required=True,
        help="Directory containing one subdirectory per ratio sample",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for PDF and TSV",
    )
    parser.add_argument("--span", type=int, default=100, help="Total plotted window size in bp")
    parser.add_argument("--min-depth", type=int, default=10, help="Minimum depth in sample and control")
    parser.add_argument(
        "--motif-pick-ratio",
        type=float,
        default=0.1,
        help="Fraction of motifs retained by score",
    )
    parser.add_argument("--pick", choices=["top", "bottom"], default="top")
    parser.add_argument("--method", choices=["median", "mean"], default="median")
    parser.add_argument(
        "--flank-exclusion",
        type=int,
        default=20,
        help="Exclude the central interval [-N, N) when calculating flank mean",
    )
    parser.add_argument(
        "--sample-regex",
        default=None,
        help=(
            "Optional custom regex used to match TF sample directory names. "
            "The default matches TF, TFFL, or TFDBD as a complete token."
        ),
    )
    parser.add_argument(
        "--sample-regex-name",
        default=None,
        help=(
            "Optional short name for the custom sample regex. It is appended to "
            "the PDF/TSV filenames and shown in the PDF title. If omitted, a "
            "filename-safe form of --sample-regex is used."
        ),
    )
    parser.add_argument(
        "--output-tag",
        default=None,
        help=(
            "Optional extra label appended to the PDF and TSV filenames. "
            "Unsafe filename characters are replaced with underscores."
        ),
    )
    parser.add_argument(
        "--raw-ylim",
        nargs=2,
        type=float,
        metavar=("LOW", "HIGH"),
        default=None,
        help="Y-axis limits for the unnormalized subplots on both pages",
    )
    parser.add_argument(
        "--normalized-ylim",
        nargs=2,
        type=float,
        metavar=("LOW", "HIGH"),
        default=None,
        help="Y-axis limits for the flank-normalized subplots on both pages",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()

    ratio_dir = args.ratio_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    if not ratio_dir.is_dir():
        raise FileNotFoundError(f"Ratio directory does not exist: {ratio_dir}")
    if args.flank_exclusion < 0 or args.flank_exclusion >= args.span // 2:
        raise ValueError("flank_exclusion must satisfy 0 <= N < span/2")

    output_dir.mkdir(parents=True, exist_ok=True)

    if args.sample_regex_name is not None and not args.sample_regex:
        raise ValueError("--sample-regex-name requires --sample-regex")
    sample_filter_label = resolve_sample_filter_label(
        sample_regex=args.sample_regex,
        sample_regex_name=args.sample_regex_name,
    )

    motif_bed = resolve_motif_bed(args.tf, args.motif_bed, args.motif_dir)
    control_dir = resolve_sample_dir(ratio_dir, args.positive_control)
    control_file = resolve_sites_file(control_dir)
    sample_files = discover_tf_samples(
        ratio_dir=ratio_dir,
        tf_name=args.tf,
        control_dir=control_dir,
        sample_regex=args.sample_regex,
    )

    print(f"[INFO] Script version: {SCRIPT_VERSION}")
    print(f"[INFO] TF: {args.tf}")
    print(f"[INFO] Ratio directory: {ratio_dir}")
    print(f"[INFO] Motif BED: {motif_bed}")
    print(f"[INFO] Positive control: {control_dir.name}")
    print(f"[INFO] Positive-control sites: {control_file}")
    if args.sample_regex:
        print(f"[INFO] Sample regex: {args.sample_regex}")
        print(f"[INFO] Sample regex output label: {sample_filter_label}")
    else:
        print("[INFO] Sample regex: not specified; using default TF matching")
    print(f"[INFO] Matched TF samples: {len(sample_files)}")
    for sample_name, sites_file in sample_files:
        print(f"       - {sample_name}\n         {sites_file}")

    motifs, score_threshold, malformed = load_and_select_motifs(
        motif_bed=motif_bed,
        motif_pick_ratio=args.motif_pick_ratio,
        pick=args.pick,
    )
    print(f"[INFO] Selected motifs: {len(motifs):,}")
    print(f"[INFO] Motif score threshold: {score_threshold:g}")
    if malformed:
        eprint(f"[WARN] Ignored malformed/invalid motif records: {malformed:,}")

    x, corrected_profiles, raw_profiles, corrected_counts, raw_counts = compute_profiles(
        sample_files=sample_files,
        control_file=control_file,
        motifs=motifs,
        span=args.span,
        min_depth=args.min_depth,
        method=args.method,
    )

    corrected_normalized_profiles: Dict[str, np.ndarray] = {}
    corrected_normalization_factors: Dict[str, float] = {}
    raw_normalized_profiles: Dict[str, np.ndarray] = {}
    raw_normalization_factors: Dict[str, float] = {}

    for sample_name in corrected_profiles:
        corrected_norm, corrected_factor = flank_normalize(
            x=x,
            profile=corrected_profiles[sample_name],
            flank_exclusion=args.flank_exclusion,
        )
        corrected_normalized_profiles[sample_name] = corrected_norm
        corrected_normalization_factors[sample_name] = corrected_factor

        raw_norm, raw_factor = flank_normalize(
            x=x,
            profile=raw_profiles[sample_name],
            flank_exclusion=args.flank_exclusion,
        )
        raw_normalized_profiles[sample_name] = raw_norm
        raw_normalization_factors[sample_name] = raw_factor

    print("[INFO] Profile summary:")
    for sample_name in corrected_profiles:
        corrected_valid_positions = int(np.isfinite(corrected_profiles[sample_name]).sum())
        raw_valid_positions = int(np.isfinite(raw_profiles[sample_name]).sum())

        corrected_factor = corrected_normalization_factors[sample_name]
        raw_factor = raw_normalization_factors[sample_name]

        corrected_factor_text = "nan" if not math.isfinite(corrected_factor) else f"{corrected_factor:.6g}"
        raw_factor_text = "nan" if not math.isfinite(raw_factor) else f"{raw_factor:.6g}"

        print(
            f"       - {sample_name}: "
            f"PC-corrected valid positions={corrected_valid_positions}/{args.span}, "
            f"PC-corrected contributions={corrected_counts[sample_name]:,}, "
            f"PC flank factor={corrected_factor_text}; "
            f"raw valid positions={raw_valid_positions}/{args.span}, "
            f"raw contributions={raw_counts[sample_name]:,}, "
            f"raw flank factor={raw_factor_text}"
        )

        if corrected_valid_positions == 0:
            eprint(
                f"[WARN] {sample_name} has no valid PC-corrected plotted positions. "
                "Check chromosome naming, tabix coordinates, depth threshold, and control coverage."
            )
        if raw_valid_positions == 0:
            eprint(
                f"[WARN] {sample_name} has no valid raw plotted positions. "
                "Check chromosome naming, tabix coordinates, and depth threshold."
            )

    prefix = sanitize_filename(args.tf)
    if sample_filter_label is not None:
        prefix = f"{prefix}__{sanitize_filename(sample_filter_label)}"
    if args.output_tag:
        prefix = f"{prefix}__{sanitize_filename(args.output_tag)}"
    output_pdf = output_dir / f"{prefix}_motif_centered_conversion.pdf"
    output_tsv = output_dir / f"{prefix}_motif_centered_conversion.tsv"

    write_profile_table(
        output_tsv=output_tsv,
        x=x,
        corrected_profiles=corrected_profiles,
        corrected_normalized_profiles=corrected_normalized_profiles,
        raw_profiles=raw_profiles,
        raw_normalized_profiles=raw_normalized_profiles,
    )
    plot_profiles(
        output_pdf=output_pdf,
        tf_name=args.tf,
        motif_bed=motif_bed,
        control_name=control_dir.name,
        sample_filter_label=sample_filter_label,
        sample_regex=args.sample_regex,
        x=x,
        corrected_profiles=corrected_profiles,
        corrected_normalized_profiles=corrected_normalized_profiles,
        raw_profiles=raw_profiles,
        raw_normalized_profiles=raw_normalized_profiles,
        motif_count=len(motifs),
        score_threshold=score_threshold,
        motif_pick_ratio=args.motif_pick_ratio,
        pick=args.pick,
        flank_exclusion=args.flank_exclusion,
        raw_ylim=parse_ylim(args.raw_ylim),
        normalized_ylim=parse_ylim(args.normalized_ylim),
    )

    print(f"[DONE] PDF: {output_pdf}")
    print(f"[DONE] TSV: {output_tsv}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        eprint("\n[ERROR] Interrupted by user.")
        raise SystemExit(130)
    except Exception as exc:
        eprint(f"[ERROR] {exc}")
        raise SystemExit(1)

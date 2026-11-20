#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  bash run_motif_centered_conversion_v4.sh TF POSITIVE_CONTROL_SAMPLE [TFBS_BED]

Arguments:
  TF                       TF name, for example CTCF, NRF1, SP1
  POSITIVE_CONTROL_SAMPLE  Complete positive-control sample directory name
  TFBS_BED                 Optional explicit BED/BED.GZ path. If omitted, the
                           script finds TF.bed case-insensitively in MOTIF_DIR.

Examples:
  bash run_motif_centered_conversion_v4.sh \
      NRF1 \
      0702-cAGM-ivtFoodie-829B4-0d1x-2min-PC_L3_Q1389W0626

  bash run_motif_centered_conversion_v4.sh \
      CTCF \
      0702-cAGM-ivtFoodie-829B4-0d1x-2min-PC_L3_Q1389W0626 \
      /dshare/xielab/analysisdata/TFC/shenk/important_data/TFmotifs/jaspar/JASPAR_split_d50/CTCF.bed

Optional environment variables:
  PYTHON                Python executable; default: python
  RATIO_DIR             Default: ../02-merge-outputs/ratio relative to this script
  MOTIF_DIR             Default JASPAR_split_d50 directory
  OUTPUT_DIR            Default: motif_centered_conversion_results under this script
  SPAN                   Default: 100
  MIN_DEPTH              Default: 10
  MOTIF_PICK_RATIO       Default: 0.1
  PICK                   top or bottom; default: top
  METHOD                 median or mean; default: median
  FLANK_EXCLUSION        Default: 20
  SAMPLE_REGEX           Optional custom regex for TF sample directory matching
  SAMPLE_REGEX_NAME      Optional short output label for SAMPLE_REGEX. If omitted,
                         a filename-safe form of SAMPLE_REGEX is used automatically.
  OUTPUT_TAG             Optional extra label appended to PDF/TSV filenames
  RAW_YLIM               Optional, for example: "0 2"
  NORMALIZED_YLIM        Optional, for example: "0.5 1.5"
EOF
}

if [[ $# -lt 2 || $# -gt 3 ]]; then
    usage >&2
    exit 2
fi

TF_NAME=$1
POSITIVE_CONTROL_SAMPLE=$2
TFBS_BED=${3:-}

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PYTHON=${PYTHON:-python}
RATIO_DIR=${RATIO_DIR:-"${SCRIPT_DIR}/../02-merge-outputs/ratio"}
MOTIF_DIR=${MOTIF_DIR:-/dshare/xielab/analysisdata/TFC/shenk/important_data/TFmotifs/jaspar/JASPAR_split_d50}
OUTPUT_DIR=${OUTPUT_DIR:-"${SCRIPT_DIR}/motif_centered_conversion_results"}
SPAN=${SPAN:-100}
MIN_DEPTH=${MIN_DEPTH:-10}
MOTIF_PICK_RATIO=${MOTIF_PICK_RATIO:-0.1}
PICK=${PICK:-top}
METHOD=${METHOD:-median}
FLANK_EXCLUSION=${FLANK_EXCLUSION:-20}

PYTHON_SCRIPT="${SCRIPT_DIR}/plot_motif_centered_conversion_v4.py"

if [[ ! -f "${PYTHON_SCRIPT}" ]]; then
    echo "[ERROR] Python script not found: ${PYTHON_SCRIPT}" >&2
    exit 1
fi
if [[ ! -d "${RATIO_DIR}" ]]; then
    echo "[ERROR] Ratio directory not found: ${RATIO_DIR}" >&2
    exit 1
fi
if ! command -v "${PYTHON}" >/dev/null 2>&1; then
    echo "[ERROR] Python executable not found: ${PYTHON}" >&2
    exit 1
fi

mkdir -p -- "${OUTPUT_DIR}"

cmd=(
    "${PYTHON}" "${PYTHON_SCRIPT}"
    --tf "${TF_NAME}"
    --positive-control "${POSITIVE_CONTROL_SAMPLE}"
    --ratio-dir "${RATIO_DIR}"
    --motif-dir "${MOTIF_DIR}"
    --output-dir "${OUTPUT_DIR}"
    --span "${SPAN}"
    --min-depth "${MIN_DEPTH}"
    --motif-pick-ratio "${MOTIF_PICK_RATIO}"
    --pick "${PICK}"
    --method "${METHOD}"
    --flank-exclusion "${FLANK_EXCLUSION}"
)

if [[ -n "${TFBS_BED}" ]]; then
    cmd+=(--motif-bed "${TFBS_BED}")
fi
if [[ -n "${SAMPLE_REGEX:-}" ]]; then
    cmd+=(--sample-regex "${SAMPLE_REGEX}")
    if [[ -n "${SAMPLE_REGEX_NAME:-}" ]]; then
        cmd+=(--sample-regex-name "${SAMPLE_REGEX_NAME}")
    fi
elif [[ -n "${SAMPLE_REGEX_NAME:-}" ]]; then
    echo '[ERROR] SAMPLE_REGEX_NAME requires SAMPLE_REGEX.' >&2
    exit 2
fi
if [[ -n "${OUTPUT_TAG:-}" ]]; then
    cmd+=(--output-tag "${OUTPUT_TAG}")
fi
if [[ -n "${RAW_YLIM:-}" ]]; then
    read -r raw_low raw_high extra <<<"${RAW_YLIM}"
    if [[ -z "${raw_low:-}" || -z "${raw_high:-}" || -n "${extra:-}" ]]; then
        echo '[ERROR] RAW_YLIM must contain exactly two numbers, for example: RAW_YLIM="0 2"' >&2
        exit 2
    fi
    cmd+=(--raw-ylim "${raw_low}" "${raw_high}")
fi
if [[ -n "${NORMALIZED_YLIM:-}" ]]; then
    read -r norm_low norm_high extra <<<"${NORMALIZED_YLIM}"
    if [[ -z "${norm_low:-}" || -z "${norm_high:-}" || -n "${extra:-}" ]]; then
        echo '[ERROR] NORMALIZED_YLIM must contain exactly two numbers, for example: NORMALIZED_YLIM="0.5 1.5"' >&2
        exit 2
    fi
    cmd+=(--normalized-ylim "${norm_low}" "${norm_high}")
fi

printf '[INFO] Command:'
printf ' %q' "${cmd[@]}"
printf '\n'

"${cmd[@]}"


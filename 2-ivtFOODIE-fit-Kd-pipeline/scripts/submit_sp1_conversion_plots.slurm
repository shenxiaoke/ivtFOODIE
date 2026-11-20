#!/usr/bin/env bash
#SBATCH --job-name=sp1_conv_plots
#SBATCH --partition=compute_new
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-4%2
#SBATCH --output=/dshare/xielab/analysisdata/TFC/shenk/analysis/WZ/BANC_FOODIE/260602_FOODIEin/20260530_WZ_FOODIEin_SP1-MLX-YY1-dbd/results/logs/slurm/sp1_conv_plots_%A_%a.out
#SBATCH --error=/dshare/xielab/analysisdata/TFC/shenk/analysis/WZ/BANC_FOODIE/260602_FOODIEin/20260530_WZ_FOODIEin_SP1-MLX-YY1-dbd/results/logs/slurm/sp1_conv_plots_%A_%a.err

set -euo pipefail

ROOT=/dshare/xielab/analysisdata/TFC/shenk/analysis/WZ/BANC_FOODIE/260602_FOODIEin
SCRIPT=$ROOT/20260530_WZ_FOODIEin_SP1-MLX-YY1-dbd/ivtFOODIE_postprocess_codeRefactoring/scripts/run_sp1_fl_linear_logo_scan.py
PY=/share/home/shenk/local/anaconda3/envs/bismark/bin/python
FL_LOGO=$ROOT/seq2FL-DBD_linearKd_model/logos/unfiltered/SP1/flank_corrected/FL_neglog10Kd/ols/linear_logo.tsv
DBD_LOGO=$ROOT/seq2FL-DBD_linearKd_model/logos/unfiltered/SP1/flank_corrected/DBD_neglog10Kd/ols/coef_logo_zero_sum.tsv
JASPAR_BED=/dshare/xielab/analysisdata/TFC/shenk/important_data/TFmotifs/jaspar/JASPAR_split_d50/SP1.bed

cd "$ROOT"

case "${SLURM_ARRAY_TASK_ID}" in
  0)
    "$PY" "$SCRIPT" \
      --tf-name SP1 \
      --results-tf SP1 \
      --logo-path "$FL_LOGO" \
      --motif-name SP1_FL_linear_logo_pos10_20 \
      --output-tag FL_neglog10Kd_ols_linear_logo_pos10_20
    ;;
  1)
    "$PY" "$SCRIPT" \
      --tf-name SP1 \
      --results-tf SP1 \
      --logo-path "$DBD_LOGO" \
      --motif-name SP1_DBD_coef_logo_pos10_20 \
      --output-tag DBD_neglog10Kd_ols_coef_logo_pos10_20
    ;;
  2)
    "$PY" "$SCRIPT" \
      --tf-name SP1-FL \
      --results-tf SP1-FL \
      --logo-path "$DBD_LOGO" \
      --motif-name SP1_DBD_coef_logo_pos10_20 \
      --output-tag DBD_neglog10Kd_ols_coef_logo_pos10_20
    ;;
  3)
    "$PY" "$SCRIPT" \
      --tf-name SP1 \
      --results-tf SP1 \
      --binding-sites-bed "$JASPAR_BED" \
      --output-tag JASPAR_SP1
    ;;
  4)
    "$PY" "$SCRIPT" \
      --tf-name SP1-FL \
      --results-tf SP1-FL \
      --binding-sites-bed "$JASPAR_BED" \
      --output-tag JASPAR_SP1
    ;;
  *)
    echo "Unknown SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
    exit 2
    ;;
esac

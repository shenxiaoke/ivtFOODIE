#!/bin/bash
set -euo pipefail

FAI_OR_CHROMSIZES=${1:? "Usage: bash run_make_10_running_Kd_spearmanr.sh genome.fa.fai|chrom.sizes"}

WINDOW_BP=${WINDOW_BP:-10}
STEP_BP=${STEP_BP:-4}
MIN_N=${MIN_N:-3}

if ! command -v bedGraphToBigWig >/dev/null 2>&1; then
    echo "[ERROR] bedGraphToBigWig not found in PATH."
    exit 1
fi

PY=./make_running_Kd_spearmanr.py

if [[ ! -s "$PY" ]]; then
    echo "[ERROR] Cannot find $PY"
    exit 1
fi

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

CHROM_SIZES="${tmpdir}/chrom.sizes"

# genome.fa.fai 是 5 列；bedGraphToBigWig 需要 2 列 chrom.sizes
awk 'BEGIN{OFS="\t"} {print $1, $2}' "$FAI_OR_CHROMSIZES" > "$CHROM_SIZES"

for fitting_dir in ../../results/*/03_fitting
do
    tf_dir=$(dirname "$fitting_dir")
    tf=$(basename "$tf_dir")

    outdir="${tf_dir}/10_running_Kd_spearmanr"
    mkdir -p "$outdir"

    echo "[INFO] Processing ${tf}"

    for mode in basic flank_corrected
    do
        input="${tf_dir}/03_fitting/${tf}.fittingKd.${mode}.tsv"

        if [[ ! -s "$input" ]]; then
            echo "[WARN] Missing file: $input"
            continue
        fi

        prefix="${tf}.${mode}"

        python "$PY" \
            --input "$input" \
            --outdir "$outdir" \
            --prefix "$prefix" \
            --chrom-sizes "$CHROM_SIZES" \
            --window-bp "$WINDOW_BP" \
            --step-bp "$STEP_BP" \
            --min-n "$MIN_N"

        for signal in minusLog10Kd spearmanr
        do
            bg="${outdir}/${prefix}.${signal}.window${WINDOW_BP}bp.step${STEP_BP}bp.center4bp.bedGraph"
            bw="${outdir}/${prefix}.${signal}.window${WINDOW_BP}bp.step${STEP_BP}bp.center4bp.bw"

            if [[ ! -s "$bg" ]]; then
                echo "[WARN] Empty bedGraph: $bg"
                continue
            fi

            sorted_bg="${bg}.sorted.tmp"

            # bedGraphToBigWig 要求 case-sensitive sorted
            # 同时去除完全重复 interval，并跳过意外重叠 interval
            LC_ALL=C sort -k1,1 -k2,2n -k3,3n "$bg" \
                | awk '
                    BEGIN{OFS="\t"}
                    NF >= 4 {
                        chr=$1
                        st=$2
                        ed=$3
                        val=$4

                        if (st < 0 || ed <= st) {
                            next
                        }

                        key = chr "\t" st "\t" ed

                        if (key in seen) {
                            next
                        }
                        seen[key] = 1

                        if (chr == last_chr && st < last_ed) {
                            next
                        }

                        print chr, st, ed, val
                        last_chr = chr
                        last_ed = ed
                    }
                ' > "$sorted_bg"

            mv "$sorted_bg" "$bg"

            bedGraphToBigWig "$bg" "$CHROM_SIZES" "$bw"

            echo "[OK] $bw"
        done
    done
done

# ----------------------------------------------------------------------
# 集中软链接所有 sample/TF 的结果目录
# ----------------------------------------------------------------------

link_root="../../results/10_running_Kd_spearmanr"
mkdir -p "$link_root"

echo "[INFO] Creating symlinks under ${link_root}/"

for src in ../../results/*/10_running_Kd_spearmanr
do
    [[ -d "$src" ]] || continue

    sample_dir=$(dirname "$src")
    sample=$(basename "$sample_dir")

    if [[ "$sample" == "10_running_Kd_spearmanr" ]]; then
        continue
    fi

    link_path="${link_root}/${sample}"
    relative_target="../${sample}/10_running_Kd_spearmanr"

    ln -sfn "$relative_target" "$link_path"

    echo "[LINK] ${link_path} -> ${relative_target}"
done

echo "[DONE] Generated bigWig files under */10_running_Kd_spearmanr/"
echo "[DONE] Symlinks created under ${link_root}/"
echo "[PARAM] WINDOW_BP=${WINDOW_BP}, STEP_BP=${STEP_BP}, MIN_N=${MIN_N}"


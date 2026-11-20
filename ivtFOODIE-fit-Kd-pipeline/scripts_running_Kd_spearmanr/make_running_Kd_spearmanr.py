#!/usr/bin/env python3

import argparse
import math
import os
from collections import defaultdict


def parse_float(x):
    try:
        v = float(x)
        if math.isfinite(v):
            return v
        return float("nan")
    except Exception:
        return float("nan")


def load_chrom_sizes(path):
    chrom_sizes = {}
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            arr = line.rstrip("\n").split("\t")
            chrom_sizes[arr[0]] = int(arr[1])
    return chrom_sizes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--window-bp", type=int, default=10)
    ap.add_argument("--step-bp", type=int, default=4)
    ap.add_argument("--min-n", type=int, default=3)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    chrom_sizes = load_chrom_sizes(args.chrom_sizes)

    # key: (chrom, window_start)
    # value:
    # [
    #   sum_minus_log10_Kd,
    #   n_valid_Kd,
    #   sum_spearmanr,
    #   n_valid_spearmanr
    # ]
    agg = defaultdict(lambda: [0.0, 0, 0.0, 0])

    n_input = 0
    n_used_kd = 0
    n_used_sp = 0
    n_skip_chr = 0
    n_skip_bad_coord = 0

    with open(args.input) as fin:
        header = fin.readline()

        for line in fin:
            if not line.strip():
                continue

            arr = line.rstrip("\n").split("\t")
            n_input += 1

            if len(arr) < 9:
                continue

            ch = arr[0]

            if ch not in chrom_sizes:
                n_skip_chr += 1
                continue

            try:
                st = int(arr[1])
                ed = int(arr[2])
            except Exception:
                n_skip_bad_coord += 1
                continue

            if st < 0 or ed <= st or st >= chrom_sizes[ch]:
                n_skip_bad_coord += 1
                continue

            # basic 和 flank_corrected 最后 6 列都是：
            # Kd  n  spearmanr  err  KdErr  deltaG
            kd = parse_float(arr[-6])
            sp = parse_float(arr[-4])

            minus_log10_kd = float("nan")
            if math.isfinite(kd) and kd > 0:
                minus_log10_kd = -math.log10(kd)

            # 只为真实有数据的位置寻找覆盖它的 10 bp sliding windows
            # 不扫描全基因组空白区域
            first_w = ((st - args.window_bp + 1) // args.step_bp) * args.step_bp
            last_w = (st // args.step_bp) * args.step_bp

            for wst in range(first_w, last_w + 1, args.step_bp):
                if wst < 0:
                    continue

                wed = wst + args.window_bp

                if wed > chrom_sizes[ch]:
                    continue

                # 当前 1 bp position 必须落入这个 10 bp window
                if not (wst <= st < wed):
                    continue

                key = (ch, wst)

                if math.isfinite(minus_log10_kd):
                    agg[key][0] += minus_log10_kd
                    agg[key][1] += 1
                    n_used_kd += 1

                if math.isfinite(sp):
                    agg[key][2] += sp
                    agg[key][3] += 1
                    n_used_sp += 1

    kd_bg = os.path.join(
        args.outdir,
        f"{args.prefix}.minusLog10Kd.window{args.window_bp}bp.step{args.step_bp}bp.center4bp.bedGraph"
    )
    sp_bg = os.path.join(
        args.outdir,
        f"{args.prefix}.spearmanr.window{args.window_bp}bp.step{args.step_bp}bp.center4bp.bedGraph"
    )

    center_len = args.step_bp
    left_trim = (args.window_bp - center_len) // 2
    right_trim = args.window_bp - center_len - left_trim

    n_out_kd = 0
    n_out_sp = 0

    with open(kd_bg, "w") as fkd, open(sp_bg, "w") as fsp:
        for ch, wst in sorted(agg.keys(), key=lambda x: (x[0], x[1])):
            wed = wst + args.window_bp

            # 10 bp window 的中心 4 bp
            out_st = wst + left_trim
            out_ed = wed - right_trim

            if out_st < 0 or out_ed > chrom_sizes[ch] or out_ed <= out_st:
                continue

            kd_sum, kd_n, sp_sum, sp_n = agg[(ch, wst)]

            if kd_n >= args.min_n:
                kd_mean = kd_sum / kd_n
                fkd.write(f"{ch}\t{out_st}\t{out_ed}\t{kd_mean:.10g}\n")
                n_out_kd += 1

            if sp_n >= args.min_n:
                sp_mean = sp_sum / sp_n
                fsp.write(f"{ch}\t{out_st}\t{out_ed}\t{sp_mean:.10g}\n")
                n_out_sp += 1

    qc = os.path.join(
        args.outdir,
        f"{args.prefix}.window{args.window_bp}bp.step{args.step_bp}bp.center4bp.qc.txt"
    )

    with open(qc, "w") as f:
        f.write(f"input\t{args.input}\n")
        f.write(f"window_bp\t{args.window_bp}\n")
        f.write(f"step_bp\t{args.step_bp}\n")
        f.write(f"center_len_bp\t{center_len}\n")
        f.write(f"min_n\t{args.min_n}\n")
        f.write(f"n_input_rows\t{n_input}\n")
        f.write(f"n_used_minus_log10_Kd_values_across_windows\t{n_used_kd}\n")
        f.write(f"n_used_spearmanr_values_across_windows\t{n_used_sp}\n")
        f.write(f"n_output_minus_log10_Kd_windows\t{n_out_kd}\n")
        f.write(f"n_output_spearmanr_windows\t{n_out_sp}\n")
        f.write(f"n_skip_chrom_not_in_chrom_sizes\t{n_skip_chr}\n")
        f.write(f"n_skip_bad_coord\t{n_skip_bad_coord}\n")


if __name__ == "__main__":
    main()

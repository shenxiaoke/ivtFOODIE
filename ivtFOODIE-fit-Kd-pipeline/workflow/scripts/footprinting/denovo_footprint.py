"""
De novo footprinting模块

功能说明：
    基于Kd阈值自动发现TF结合足迹（footprints）

生物学原理：
    - TF结合位点通常包含多个连续的高亲和力碱基
    - 使用双层阈值策略：
      * 种子区域（Seed）：高置信度的结合核心（如top 5%）
      * 扩展区域（Extended）：从种子向两侧扩展（如top 10%）
    - 自动合并重叠的footprints

算法原理：
    1. 计算Kd的-log10值（Kd越小，-log10(Kd)越大）
    2. 确定种子阈值（如95th percentile）
    3. 从每个种子位点向两侧扩展，纳入超过扩展阈值的碱基
    4. 当连续遇到低于阈值的碱基时停止扩展
    5. 只保留长度在[min_len, max_len]范围内的footprints
    6. 合并重叠的footprints

参考：
    原始实现：footprinting.py
    作者：沈柯
"""

from typing import Optional, Tuple, List
from pathlib import Path
import pandas as pd
import numpy as np
from numpy.typing import NDArray
import multiprocessing
from functools import partial


def ratio_split(ratio_df: pd.DataFrame) -> dict:
    """
    按染色体分割ratio DataFrame

    加速查询：将DataFrame按染色体分组，避免每次都全表扫描

    参数：
        ratio_df: Kd拟合结果DataFrame

    返回：
        {chromosome: sub_dataframe} 字典
    """
    ratio_dict = {}
    for ch in ratio_df['ch'].unique():
        ratio_dict[ch] = ratio_df[ratio_df['ch'] == ch]
    return ratio_dict


def denovo_footprinting_single_seed(
    index: int,
    ratio_df: pd.DataFrame,
    ratio_df_dict: dict,
    kd_log_thres_seed: float,
    kd_log_thres_extend: float,
    min_len: int = 10,
    max_len: int = 30,
    max_gap: int = 5,
    tolerance: float = 0.1,
    kd_column: str = 'Kd_log'
) -> Optional[List]:
    """
    从单个种子位点扩展footprint

    算法：
        1. 从种子位点开始，向左右两侧扩展
        2. 每次扩展检查：
           - 间隙是否超过max_gap
           - 扩展后低于阈值的位点比例是否超过tolerance
        3. 如果两侧都无法扩展，停止
        4. 返回扩展后的footprint区间

    参数：
        index: 种子位点在ratio_df中的索引
        ratio_df: 完整的Kd拟合结果DataFrame
        ratio_df_dict: 按染色体分组的DataFrame字典
        kd_log_thres_seed: 种子区域的Kd阈值（-log10(Kd)）
        kd_log_thres_extend: 扩展区域的Kd阈值（-log10(Kd)）
        min_len: 最小footprint长度 (bp)
        max_len: 最大footprint长度 (bp)
        max_gap: 允许的最大间隙 (bp)
        tolerance: 允许的低于阈值位点比例（默认10%）
        kd_column: Kd列名

    返回：
        [ch, st, ed] 或 None（如果不满足长度要求）
    """
    # 提取种子位点信息
    ch = ratio_df.iloc[index]['ch']
    st = int(ratio_df.iloc[index]['st'])

    # 获取该染色体的数据，并限制在种子位点±max_len范围内
    selected_ratio_df = ratio_df_dict[ch]
    selected_ratio_df = selected_ratio_df[
        (selected_ratio_df['st'] >= st - max_len) &
        (selected_ratio_df['ed'] < st + max_len)
    ]

    # 提取位置和Kd_log值
    st_list = list(selected_ratio_df['st'])
    kd_log_list = np.array(selected_ratio_df[kd_column])

    # 找到种子位点在st_list中的索引
    try:
        center_i = st_list.index(st)
    except ValueError:
        # 种子位点不在选定区域内（不应该发生）
        return None

    # 初始化扩展标志
    stop_left, stop_right = False, False
    dormant_left, dormant_right = False, False
    left_i, right_i = center_i, center_i

    # 双向扩展
    while not stop_left or not stop_right:
        # 向左扩展
        if not stop_left:
            # 检查是否到达边界
            if left_i - 1 < 0:
                stop_left = True
            # 检查间隙是否过大
            elif st_list[left_i] - st_list[left_i - 1] > max_gap:
                stop_left = True
            # 检查扩展后低于阈值的位点比例
            elif (
                len(kd_log_list[left_i - 1:right_i + 1][
                    kd_log_list[left_i - 1:right_i + 1] < kd_log_thres_extend
                ]) > tolerance * (right_i - left_i + 2)
            ):
                # 太多位点低于阈值，标记为休眠
                dormant_left = True
                if dormant_right or stop_right:
                    break  # 两侧都无法扩展，停止
            else:
                # 可以扩展，重置休眠状态
                dormant_left, dormant_right = False, False
                left_i -= 1

        # 向右扩展
        if not stop_right:
            # 检查是否到达边界
            if right_i + 1 >= len(st_list):
                stop_right = True
            # 检查间隙是否过大
            elif st_list[right_i + 1] - st_list[right_i] > max_gap:
                stop_right = True
            # 检查扩展后低于阈值的位点比例
            elif (
                len(kd_log_list[left_i:right_i + 2][
                    kd_log_list[left_i:right_i + 2] < kd_log_thres_extend
                ]) > tolerance * (right_i - left_i + 2)
            ):
                # 太多位点低于阈值，标记为休眠
                dormant_right = True
                if dormant_left or stop_left:
                    break  # 两侧都无法扩展，停止
            else:
                # 可以扩展，重置休眠状态
                dormant_left, dormant_right = False, False
                right_i += 1

    # 检查footprint长度是否满足要求
    footprint_length = st_list[right_i] + 1 - st_list[left_i]
    if min_len <= footprint_length <= max_len:
        return [ch, st_list[left_i], st_list[right_i] + 1]

    return None


def process_rows_parallel(
    ratio_df: pd.DataFrame,
    kd_log_thres_seed: float,
    kd_log_thres_extend: float,
    num_processes: Optional[int] = None,
    kd_column: str = 'Kd_log',
    min_len: int = 10,
    max_len: int = 30,
    max_gap: int = 5,
    tolerance: float = 0.1
) -> pd.DataFrame:
    """
    并行处理所有种子位点，生成footprints

    流程：
        1. 筛选种子位点（Kd_log >= kd_log_thres_seed）
        2. 对每个种子位点并行扩展footprint
        3. 合并重叠的footprints
        4. 返回footprint DataFrame

    参数：
        ratio_df: Kd拟合结果DataFrame（必须包含'ch', 'st', 'ed', Kd_column列）
        kd_log_thres_seed: 种子区域阈值（-log10(Kd)）
        kd_log_thres_extend: 扩展区域阈值（-log10(Kd)）
        num_processes: 并行进程数（None表示使用所有CPU核心）
        kd_column: Kd列名
        min_len: 最小footprint长度 (bp)
        max_len: 最大footprint长度 (bp)
        max_gap: 允许的最大间隙 (bp)
        tolerance: 允许的低于阈值位点比例

    返回：
        Footprint DataFrame (columns: ch, st, ed)
    """
    # 按染色体分组（加速查询）
    ratio_df_dict = ratio_split(ratio_df)

    # 筛选种子位点
    seed_indices = list(ratio_df[ratio_df[kd_column] >= kd_log_thres_seed].index)
    print(f'Seed sites: {len(seed_indices)}')

    # 创建partial函数（固定参数）
    partial_footprinting = partial(
        denovo_footprinting_single_seed,
        ratio_df=ratio_df,
        ratio_df_dict=ratio_df_dict,
        kd_log_thres_seed=kd_log_thres_seed,
        kd_log_thres_extend=kd_log_thres_extend,
        min_len=min_len,
        max_len=max_len,
        max_gap=max_gap,
        tolerance=tolerance,
        kd_column=kd_column
    )

    # 确定并行进程数
    if num_processes is None:
        num_processes = multiprocessing.cpu_count()

    print(f'Using {num_processes} processes')

    # 并行处理所有种子位点
    with multiprocessing.Pool(processes=num_processes) as pool:
        footprints = pool.map(partial_footprinting, seed_indices)

    # 过滤None结果
    footprints = [fp for fp in footprints if fp is not None]
    print(f'Raw footprints: {len(footprints)}')

    # 按染色体和位置排序
    footprints = sorted(footprints, key=lambda t: (t[0], t[1]))

    # 合并重叠的footprints
    ch_list, st_list, ed_list = [], [], []
    prev_ch, prev_st, prev_ed = ('chr1', 0, 0)

    for footprint in footprints:
        ch, st, ed = footprint

        # 检查是否与前一个footprint重叠
        if ch == prev_ch and st <= prev_ed:
            # 重叠，跳过（保留前一个）
            continue

        # 不重叠，加入结果
        prev_ch, prev_st, prev_ed = ch, st, ed
        ch_list.append(ch)
        st_list.append(st)
        ed_list.append(ed)

    # 构建DataFrame
    footprint_df = pd.DataFrame({
        'ch': ch_list,
        'st': st_list,
        'ed': ed_list
    })

    # 确保数据类型正确
    footprint_df['st'] = footprint_df['st'].astype(np.int32)
    footprint_df['ed'] = footprint_df['ed'].astype(np.int32)

    print(f'Merged footprints: {len(footprint_df)}')

    return footprint_df


def denovo_footprinting_from_kd_quantile(
    kd_df: pd.DataFrame,
    seed_quantile: float = 0.95,
    extend_quantile: float = 0.90,
    num_processes: Optional[int] = None,
    min_len: int = 10,
    max_len: int = 30,
    max_gap: int = 5,
    tolerance: float = 0.1,
    kd_column: str = 'Kd'
) -> pd.DataFrame:
    """
    从Kd分位数自动确定阈值并生成footprints

    这是Snakemake workflow中使用的主要接口函数

    参数：
        kd_df: Kd拟合结果DataFrame
        seed_quantile: 种子区域Kd分位数（默认0.95 = top 5%）
        extend_quantile: 扩展区域Kd分位数（默认0.90 = top 10%）
        num_processes: 并行进程数
        min_len: 最小footprint长度 (bp)
        max_len: 最大footprint长度 (bp)
        max_gap: 允许的最大间隙 (bp)
        tolerance: 允许的低于阈值位点比例
        kd_column: Kd列名

    返回：
        Footprint DataFrame
    """
    # 计算Kd_log（-log10(Kd)）
    kd_df['Kd_log'] = -np.log10(kd_df[kd_column])

    # 计算阈值
    kd_log_thres_seed = kd_df['Kd_log'].quantile(seed_quantile)
    kd_log_thres_extend = kd_df['Kd_log'].quantile(extend_quantile)

    print(f"\n=== De novo Footprinting ===")
    print(f"Seed quantile: {seed_quantile} (top {100*(1-seed_quantile):.1f}%)")
    print(f"Extend quantile: {extend_quantile} (top {100*(1-extend_quantile):.1f}%)")
    print(f"Seed threshold: Kd_log >= {kd_log_thres_seed:.2f} (Kd <= {10**(-kd_log_thres_seed):.2e} M)")
    print(f"Extend threshold: Kd_log >= {kd_log_thres_extend:.2f} (Kd <= {10**(-kd_log_thres_extend):.2e} M)")
    print(f"Footprint length: [{min_len}, {max_len}] bp")
    print(f"Max gap: {max_gap} bp")
    print(f"Tolerance: {tolerance*100:.0f}%\n")

    # 运行footprinting
    footprint_df = process_rows_parallel(
        ratio_df=kd_df,
        kd_log_thres_seed=kd_log_thres_seed,
        kd_log_thres_extend=kd_log_thres_extend,
        num_processes=num_processes,
        kd_column='Kd_log',
        min_len=min_len,
        max_len=max_len,
        max_gap=max_gap,
        tolerance=tolerance
    )

    return footprint_df


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='De novo footprinting: 从Kd值自动发现TF结合足迹'
    )

    parser.add_argument('--in_f', type=str, required=True, help='输入Kd文件（TSV格式）')
    parser.add_argument('--out_f', type=str, required=True, help='输出footprint文件（BED格式）')
    parser.add_argument('--seed_quantile', type=float, default=0.95,
                        help='种子区域Kd分位数（默认0.95 = top 5%%）')
    parser.add_argument('--extend_quantile', type=float, default=0.90,
                        help='扩展区域Kd分位数（默认0.90 = top 10%%）')
    parser.add_argument('--num_processes', type=int, default=None,
                        help='并行进程数（默认使用所有CPU核心）')
    parser.add_argument('--min_len', type=int, default=10,
                        help='最小footprint长度 (bp)')
    parser.add_argument('--max_len', type=int, default=30,
                        help='最大footprint长度 (bp)')
    parser.add_argument('--max_gap', type=int, default=5,
                        help='允许的最大间隙 (bp)')
    parser.add_argument('--tolerance', type=float, default=0.1,
                        help='允许的低于阈值位点比例（默认0.1 = 10%%）')

    args = parser.parse_args()

    print(f"Input: {args.in_f}")
    print(f"Output: {args.out_f}")

    # 读取Kd数据
    print("\nLoading Kd data...")
    kd_df = pd.read_csv(args.in_f, sep='\t')
    print(f"Loaded {len(kd_df)} sites")

    # Footprinting
    footprint_df = denovo_footprinting_from_kd_quantile(
        kd_df=kd_df,
        seed_quantile=args.seed_quantile,
        extend_quantile=args.extend_quantile,
        num_processes=args.num_processes,
        min_len=args.min_len,
        max_len=args.max_len,
        max_gap=args.max_gap,
        tolerance=args.tolerance
    )

    # 输出BED文件（无表头）
    Path(args.out_f).parent.mkdir(parents=True, exist_ok=True)
    footprint_df.to_csv(args.out_f, sep='\t', header=False, index=False)

    print(f"\n=== Complete ===")
    print(f"Output saved to: {args.out_f}")
    print(f"Total footprints: {len(footprint_df)}")

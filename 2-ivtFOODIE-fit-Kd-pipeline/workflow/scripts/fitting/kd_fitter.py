"""
Kd拟合模块（支持Flank校正）

功能说明：
    从conversion ratio数据拟合TF-DNA结合的解离常数 (Kd)

生物学原理：
    - 使用Hill方程描述TF浓度与DNA结合的关系
    - Kd越小，表示亲和力越强
    - 可选Flank校正：消除DNA序列上下文对conversion的影响

数学模型：
    binding_fraction = 1 / (1 + (Kd/[TF])^n)
    conversion_ratio = 1 - binding_fraction

Flank校正原理：
    - 即使没有TF结合，不同序列上下文的C碱基对APOBEC3A的敏感性不同
    - 通过已知的高亲和力motif，计算其flanking区域的平均conversion ratio
    - 用flanking ratio作为baseline，校正每个样本的数据
"""

from typing import List, Optional, Tuple
from pathlib import Path
import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy import stats
from scipy.optimize import curve_fit

# 导入已有的Hill方程模块
import sys
sys.path.append(str(Path(__file__).parent.parent / "fitting"))
from hill_equation import hill_equation, hill_equation_n1


def fit_kd_single_site(
    ratios: NDArray[np.float64],
    concentrations: NDArray[np.float64],
    pc_ratio: float,
    fix_n_to_1: bool = True
) -> Tuple[float, float, float, float, float]:
    """
    单个C位点的Kd拟合

    参数：
        ratios: 各浓度下的conversion ratio数组
        concentrations: TF浓度数组 (M)
        pc_ratio: PC对照的conversion ratio（用于归一化）
        fix_n_to_1: 是否固定Hill系数为1

    返回：
        (Kd, n, spearman_correlation, error, KdErr)
            Kd: 解离常数 (M)
            n: Hill系数（若fix_n_to_1=True则固定为1）
            spearman_correlation: Spearman相关系数（浓度与ratio的负相关）
            error: 拟合误差（MSE * 100）
            KdErr: Kd的相对误差
    """
    # 计算Spearman相关系数（浓度越高，conversion ratio应该越低）
    spearman_corr = -stats.spearmanr(concentrations, ratios)[0]

    # PC校正
    if pc_ratio == 0:
        return np.nan, np.nan, spearman_corr, np.nan, np.nan

    ratios_normalized = ratios / pc_ratio

    # 转换为binding fraction: binding = 1 - conversion
    binding_fractions = 1.0 - ratios_normalized

    # 裁剪到[0, 1]范围（避免数值问题）
    binding_fractions = np.clip(binding_fractions, 0, 1)

    # 拟合Hill方程
    try:
        if fix_n_to_1:
            # 固定n=1，只拟合Kd
            popt, pcov = curve_fit(
                lambda c, kd: hill_equation_n1(c, kd),
                concentrations,
                binding_fractions,
                p0=[1e-7]  # 初始Kd猜测值 100nM
            )
            kd = popt[0]
            n = 1.0

            # 计算预测值和误差
            predicted = hill_equation_n1(concentrations, kd)

        else:
            # 同时拟合Kd和n
            popt, pcov = curve_fit(
                lambda c, kd, n: hill_equation(c, kd, n),
                concentrations,
                binding_fractions,
                p0=[1e-7, 1],  # 初始Kd=100nM, n=1
                bounds=((0, 0), (1, 10))  # Kd∈(0,1M), n∈(0,10)
            )
            kd, n = popt

            # 计算预测值和误差
            predicted = hill_equation(concentrations, kd, n)

        # 计算拟合误差（MSE）
        mse = np.mean((predicted - binding_fractions) ** 2)
        error = 100 * mse

        # 计算Kd的相对误差
        kd_err = np.sqrt(np.diag(pcov))[0] / kd

        return kd, n, spearman_corr, error, kd_err

    except RuntimeError:
        # 拟合失败
        return np.nan, 1.0, spearman_corr, np.nan, np.nan


class KdFitter:
    """
    全基因组Kd拟合器（支持Flank校正）
    """

    def __init__(
        self,
        ratio_df: pd.DataFrame,
        concentrations: List[float],
        fitting_indices: List[int],
        fix_n_to_1: bool = True,
        pc_column: str = 'ratio_adjust'
    ):
        """
        初始化KdFitter

        参数：
            ratio_df: 样本ratio数据表（来自integration步骤）
            concentrations: TF浓度序列 (M)
            fitting_indices: 用于拟合的样本索引（0-based）
            fix_n_to_1: 是否固定Hill系数为1
            pc_column: PC对照的列名
        """
        self.ratio_df = ratio_df.copy()
        self.concentrations = np.array(concentrations)
        self.fitting_indices = np.array(fitting_indices)
        self.fix_n_to_1 = fix_n_to_1
        self.pc_column = pc_column

        # 添加计数器列（用于打印进度）
        self.ratio_df['counter'] = self.ratio_df.index

    def _compute_flank_conversion(
        self,
        motif_bed: Path,
        motif_pick_ratio: float = 0.1,
        span: int = 100
    ) -> List[NDArray[np.float64]]:
        """
        计算motif区域的flank conversion profile

        原理：
            1. 选择高分motif（top motif_pick_ratio）
            2. 对每个motif，提取其中心±span范围内的conversion ratio
            3. 计算每个相对位置的中位数conversion ratio
            4. 返回各样本的flank conversion profile

        参数：
            motif_bed: Motif位点BED文件路径
            motif_pick_ratio: 使用top多少比例的motif（默认0.1 = top 10%）
            span: Motif中心±span的范围 (bp)

        返回：
            List[Array]: 每个样本的flank conversion profile (长度 = 2*span)
        """
        # 按染色体分组ratio_df（加速查询）
        ratio_dict = {}
        for ch in self.ratio_df['ch'].unique():
            ratio_dict[ch] = self.ratio_df[self.ratio_df['ch'] == ch]

        # 读取motif文件并筛选高分motif
        motifs = []
        motif_scores = []
        has_score_column = False

        with open(motif_bed, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                motifs.append(parts)

                # 检查是否有第5列（score列）
                if len(parts) >= 5:
                    try:
                        score = float(parts[4])
                        motif_scores.append(score)
                        has_score_column = True
                    except ValueError:
                        # 第5列不是数字，当作没有score列处理
                        has_score_column = False
                        break

        # 根据是否有score列选择motifs
        if has_score_column and len(motif_scores) > 0:
            # 有score列：选择top motifs
            score_threshold = np.percentile(motif_scores, 100 * (1 - motif_pick_ratio))
            chosen_motifs = [m for m, s in zip(motifs, motif_scores) if s >= score_threshold]
            print(f"Using {len(chosen_motifs)} / {len(motifs)} motifs with scores (top {motif_pick_ratio*100}%)")
            print(f"Score threshold: {score_threshold:.4f}")
        else:
            # 没有score列：使用所有motifs或按位置随机采样
            if motif_pick_ratio < 1.0:
                # 随机采样指定比例的motifs
                n_chosen = max(1, int(len(motifs) * motif_pick_ratio))
                chosen_indices = np.random.choice(len(motifs), n_chosen, replace=False)
                chosen_motifs = [motifs[i] for i in chosen_indices]
                print(f"Warning: Motif BED file has no score column (only {len(motifs[0])} columns)")
                print(f"Randomly sampling {len(chosen_motifs)} / {len(motifs)} motifs ({motif_pick_ratio*100}%)")
            else:
                # 使用全部motifs
                chosen_motifs = motifs
                print(f"Warning: Motif BED file has no score column (only {len(motifs[0])} columns)")
                print(f"Using all {len(chosen_motifs)} motifs")

        # 找出所有样本列
        sample_columns = [
            col for col in self.ratio_df.columns
            if col.startswith('ratio_') and col != self.pc_column
        ]

        # 初始化conversion lists
        conv_lists = [[[] for _ in range(2 * span)] for __ in range(len(sample_columns))]

        # 遍历每个motif，收集flank conversion
        for motif_idx, motif in enumerate(chosen_motifs):
            if motif_idx % 100 == 0:
                print(f"Processing motif {motif_idx}/{len(chosen_motifs)}", end='\r')

            ch = motif[0]
            st = int(motif[1])
            ed = int(motif[2])
            center = (st + ed) // 2

            # 检查该染色体是否有数据
            if ch not in ratio_dict:
                continue

            # 查询该motif附近的C位点
            selected_df = ratio_dict[ch]
            selected_df = selected_df[
                (selected_df['st'] >= st - span) &
                (selected_df['st'] < ed + span)
            ]

            # 遍历附近的每个C位点
            for _, row in selected_df.iterrows():
                c_pos = int(row['st'])
                relative_pos = c_pos - center + span  # 相对于motif中心的位置

                # 检查PC ratio
                pc_ratio = row[self.pc_column]
                if pc_ratio == 0:
                    continue

                # 只保留在范围内的位点
                if 0 <= relative_pos < 2 * span:
                    # 对每个样本，记录conversion/PC_ratio
                    for sample_idx, sample_col in enumerate(sample_columns):
                        conversion = row[sample_col]
                        conv_lists[sample_idx][relative_pos].append(conversion / pc_ratio)

        print()  # 换行

        # 计算每个位置的中位数
        flank_profiles = []
        for sample_conv_list in conv_lists:
            profile = np.array([
                np.nanmedian(position_values) if len(position_values) > 0 else np.nan
                for position_values in sample_conv_list
            ])
            flank_profiles.append(profile)

        return flank_profiles

    def _apply_flank_correction(
        self,
        flank_profiles: List[NDArray[np.float64]],
        flank_size: int = 50
    ):
        """
        应用flank校正到ratio_df

        原理：
            - 使用motif两侧的flank区域（motif center ± [flank_size, span]）的平均conversion作为baseline
            - 将每个样本的ratio除以对应的flank baseline

        参数：
            flank_profiles: 各样本的flank conversion profile
            flank_size: Flank区域大小（只使用profile两端各flank_size bp作为baseline）
        """
        # 找出所有样本列
        sample_columns = [
            col for col in self.ratio_df.columns
            if col.startswith('ratio_') and col != self.pc_column
        ]

        # 计算每个样本的flank baseline（两端各flank_size bp的平均值）
        flank_baselines = []
        for profile in flank_profiles:
            # 两端各取flank_size bp
            left_flank = profile[:flank_size]
            right_flank = profile[-flank_size:]
            baseline = (np.nanmean(left_flank) + np.nanmean(right_flank)) / 2.0
            flank_baselines.append(baseline)

        print("Flank baselines:")
        for i, baseline in enumerate(flank_baselines):
            print(f"  Sample {i+1}: {baseline:.4f}")

        # 对每个样本应用校正
        for sample_idx, sample_col in enumerate(sample_columns):
            baseline = flank_baselines[sample_idx]
            corrected_col = f"{sample_col}_flank_adjust"
            self.ratio_df[corrected_col] = self.ratio_df[sample_col] / baseline

    def fit_all_sites(
        self,
        motif_bed: Optional[Path] = None,
        motif_pick_ratio: float = 0.1,
        flank_size: int = 50,
        motif_span: int = 100
    ) -> pd.DataFrame:
        """
        拟合所有C位点的Kd

        参数：
            motif_bed: Motif BED文件路径（如果提供，则进行flank校正）
            motif_pick_ratio: 使用top多少比例的motif进行flank校正
            flank_size: Flank区域大小（用于计算baseline）
            motif_span: Motif中心±span的范围

        返回：
            添加了Kd拟合结果的DataFrame
        """
        # 如果提供了motif文件，进行flank校正
        use_flank_correction = motif_bed is not None

        if use_flank_correction:
            print("\n=== Flank Correction ===")
            flank_profiles = self._compute_flank_conversion(
                motif_bed=motif_bed,
                motif_pick_ratio=motif_pick_ratio,
                span=motif_span
            )

            # 添加小量避免除零
            flank_profiles = [profile + 1e-6 for profile in flank_profiles]

            # 应用校正
            self._apply_flank_correction(
                flank_profiles=flank_profiles,
                flank_size=flank_size
            )

            print("Flank correction applied.\n")

        # 移除NA
        self.ratio_df = self.ratio_df.dropna().copy()

        print(f"\n=== Kd Fitting ===")
        print(f"Total sites to fit: {len(self.ratio_df)}")
        print(f"Fix n to 1: {self.fix_n_to_1}")
        print(f"Concentrations: {self.concentrations}")
        print(f"Fitting indices: {self.fitting_indices}\n")

        # 准备拟合
        def fit_row(row):
            """对单行数据进行Kd拟合"""
            if row['counter'] % 1000 == 0:
                print(f"Processed {row['counter']} / {len(self.ratio_df)} sites", end='\r')

            # 提取ratio数据
            if use_flank_correction:
                # 使用校正后的数据
                ratios = np.array([
                    row[f'ratio_{i+1}_flank_adjust']
                    for i in self.fitting_indices
                ])
            else:
                # 使用原始数据
                ratios = np.array([
                    row[f'ratio_{i+1}']
                    for i in self.fitting_indices
                ])

            # 提取对应的浓度
            used_concentrations = self.concentrations[self.fitting_indices]

            # PC ratio
            pc_ratio = row[self.pc_column]

            # 拟合
            result = fit_kd_single_site(
                ratios=ratios,
                concentrations=used_concentrations,
                pc_ratio=pc_ratio,
                fix_n_to_1=self.fix_n_to_1
            )

            return result

        # 执行拟合
        self.ratio_df['Kd_n_spearmanr_err_KdErr'] = self.ratio_df.apply(
            fit_row, axis=1
        )

        print()  # 换行

        # 拆分结果到独立列
        self.ratio_df['Kd'] = self.ratio_df['Kd_n_spearmanr_err_KdErr'].apply(lambda x: x[0])
        self.ratio_df['n'] = self.ratio_df['Kd_n_spearmanr_err_KdErr'].apply(lambda x: x[1])
        self.ratio_df['spearmanr'] = self.ratio_df['Kd_n_spearmanr_err_KdErr'].apply(lambda x: x[2])
        self.ratio_df['err'] = self.ratio_df['Kd_n_spearmanr_err_KdErr'].apply(lambda x: x[3])
        self.ratio_df['KdErr'] = self.ratio_df['Kd_n_spearmanr_err_KdErr'].apply(lambda x: x[4])

        # 删除临时列
        self.ratio_df.drop('Kd_n_spearmanr_err_KdErr', axis=1, inplace=True)

        # 计算ΔG (Gibbs free energy)
        R = 8.31446261815324  # J·K⁻¹·mol⁻¹
        T = 300  # K
        self.ratio_df['deltaG'] = -R * T * np.log(self.ratio_df['Kd']) / 1e3  # kJ/mol

        print("\n=== Fitting Complete ===")
        valid_kd_count = self.ratio_df['Kd'].notna().sum()
        print(f"Valid Kd fits: {valid_kd_count} / {len(self.ratio_df)}")
        print(f"Success rate: {100 * valid_kd_count / len(self.ratio_df):.1f}%")

        return self.ratio_df


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='拟合TF-DNA结合的Kd值')

    parser.add_argument('--in_f', type=str, required=True, help='输入ratio文件')
    parser.add_argument('--out_f', type=str, required=True, help='输出Kd文件')
    parser.add_argument('--concentration', type=str, required=True,
                        help='TF浓度序列（逗号分隔），单位M')
    parser.add_argument('--fitting_used_index', type=str, required=True,
                        help='用于拟合的样本索引（逗号分隔，0-based）')
    parser.add_argument('--fix_n_to_1', type=str, default='True',
                        help='是否固定Hill系数为1 (True/False)')
    parser.add_argument('--motif_f', type=str, default=None,
                        help='Motif BED文件路径（如果提供，则进行flank校正）')
    parser.add_argument('--PC_name', type=str, default='ratio_adjust',
                        help='PC对照的列名')
    parser.add_argument('--motif_pick_ratio', type=float, default=0.1,
                        help='Flank校正使用的top motif比例')
    parser.add_argument('--flank_size', type=int, default=50,
                        help='Flank区域大小 (bp)')
    parser.add_argument('--motif_span', type=int, default=100,
                        help='Motif中心±span的范围 (bp)')

    args = parser.parse_args()

    # 解析参数
    concentrations = list(map(float, args.concentration.split(',')))
    fitting_indices = list(map(int, args.fitting_used_index.split(',')))
    fix_n_to_1 = args.fix_n_to_1.lower() in ['true', '1', 'yes']

    print(f"Input: {args.in_f}")
    print(f"Output: {args.out_f}")
    print(f"Concentrations: {concentrations}")
    print(f"Fitting indices: {fitting_indices}")
    print(f"Fix n to 1: {fix_n_to_1}")
    print(f"Motif file: {args.motif_f}")
    print(f"PC column: {args.PC_name}")

    # 读取数据
    print("\nLoading ratio data...")
    ratio_df = pd.read_csv(args.in_f, sep='\t')
    print(f"Loaded {len(ratio_df)} sites")

    # 创建fitter
    fitter = KdFitter(
        ratio_df=ratio_df,
        concentrations=concentrations,
        fitting_indices=fitting_indices,
        fix_n_to_1=fix_n_to_1,
        pc_column=args.PC_name
    )

    # 拟合
    result_df = fitter.fit_all_sites(
        motif_bed=Path(args.motif_f) if args.motif_f else None,
        motif_pick_ratio=args.motif_pick_ratio,
        flank_size=args.flank_size,
        motif_span=args.motif_span
    )

    # 输出
    Path(args.out_f).parent.mkdir(parents=True, exist_ok=True)
    result_df.dropna(subset=['Kd']).to_csv(args.out_f, sep='\t', index=False, na_rep='NULL')

    print(f"\nResults saved to: {args.out_f}")

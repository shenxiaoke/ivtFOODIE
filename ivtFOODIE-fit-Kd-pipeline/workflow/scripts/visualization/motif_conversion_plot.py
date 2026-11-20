"""
Motif-centered conversion rate plot script

功能：
    基于motifCenteredPlot_siteWise类实现，
    绘制不同浓度下的转化率图
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # 非交互式后端
import matplotlib.pyplot as plt
import pysam


class MotifCenteredPlotSiteWise:
    """
    Motif中心化转化率绘图类
    
    参考：old/functions.py:201
    """
    
    def __init__(self, gz_f, gz_control_f, motif_f, motif_pick_ratio=0.1, span=100, min_depth=10):
        """
        初始化
        
        参数：
            gz_f: 样本ratio文件（tabix压缩）
            gz_control_f: PC对照文件（tabix压缩）
            motif_f: Motif结合位点BED文件
            motif_pick_ratio: 选择top motif的比例
            span: Motif中心化窗口大小 (bp)
            min_depth: 最小测序深度
        """
        self.gz = pysam.TabixFile(gz_f)
        self.gz_control = pysam.TabixFile(gz_control_f)
        self.motif_f = motif_f
        self.motif_pick_ratio = motif_pick_ratio
        self.span = span
        self.min_depth = min_depth
        self.chosen_motifs = None
        self.conversion = None
    
    def motifs_list(self, pick_top_or_bottom='top', motif_pick_ratio=None):
        """
        选择motif列表
        
        参数：
            pick_top_or_bottom: 'top' 或 'bottom'
            motif_pick_ratio: 选择比例（覆盖初始化时的值）
        
        返回：
            选择的motif列表
        """
        if motif_pick_ratio is not None:
            self.motif_pick_ratio = motif_pick_ratio
        
        motifs = []
        motif_score_list = []
        count = 0
        
        with open(self.motif_f, 'r') as f:
            for l in f:
                count += 1
                if count % 1000 == 0:
                    print(f"Reading motifs: {count}", end='\r', file=sys.stderr)
                
                l = l.strip().split('\t')
                if len(l) < 5:
                    continue
                
                motifs.append(l)
                try:
                    motif_score_list.append(float(l[4]))
                except (ValueError, IndexError):
                    motif_score_list.append(0.0)
        
        print(f"\nTotal motifs read: {len(motifs)}", file=sys.stderr)
        
        if pick_top_or_bottom == 'top':
            self.score_thres = np.percentile(motif_score_list, 100 * (1 - self.motif_pick_ratio))
            self.chosen_motifs = [
                m for m in motifs 
                if len(m) >= 5 and float(m[4]) >= self.score_thres
            ]
        elif pick_top_or_bottom == 'bottom':
            self.score_thres = np.percentile(motif_score_list, 100 * self.motif_pick_ratio)
            self.chosen_motifs = [
                m for m in motifs 
                if len(m) >= 5 and float(m[4]) <= self.score_thres
            ]
        else:
            print('pick_top_or_bottom para error', file=sys.stderr)
            self.chosen_motifs = []
        
        print(f'Chosen motifs number: {len(self.chosen_motifs)}', file=sys.stderr)
        print(f'Score threshold: {self.score_thres}', file=sys.stderr)
        
        return self.chosen_motifs
    
    def conversion_profile(self, method='median'):
        """
        计算转化率profile
        
        参数：
            method: 'median' 或 'mean'
        
        返回：
            转化率数组
        """
        if self.chosen_motifs is None or len(self.chosen_motifs) == 0:
            print("Warning: No motifs selected. Returning empty profile.", file=sys.stderr)
            self.conversion = np.full(self.span, np.nan)
            return self.conversion
        
        self.conversion_list = [[] for _ in range(self.span)]
        count = 0
        
        for motif in self.chosen_motifs:
            count += 1
            if count % 1000 == 0:
                print(f"Processing motifs: {count}/{len(self.chosen_motifs)}", end='\r', file=sys.stderr)
            
            if len(motif) < 6:
                continue
            
            ch, st, ed = motif[0], int(motif[1]), int(motif[2])
            motif_name = motif[3] if len(motif) > 3 else 'motif'
            score = motif[4] if len(motif) > 4 else '0'
            strand = motif[5] if len(motif) > 5 else '+'
            
            center = (st + ed) // 2
            
            try:
                records = list(self.gz.fetch(ch, max(0, center - self.span // 2), center + self.span // 2))
                records_control = list(self.gz_control.fetch(ch, max(0, center - self.span // 2), center + self.span // 2))
            except Exception as e:
                continue
            
            # 解析records
            records_dict = {}
            for record in records:
                fields = record.split('\t')
                if len(fields) >= 6:
                    pos = int(fields[1])
                    records_dict[pos] = fields
            
            records_control_dict = {}
            for record in records_control:
                fields = record.split('\t')
                if len(fields) >= 6:
                    pos = int(fields[1])
                    records_control_dict[pos] = fields
            
            # 计算转化率
            for pos in set(records_dict.keys()) & set(records_control_dict.keys()):
                record = records_dict[pos]
                record_control = records_control_dict[pos]
                
                try:
                    depth = int(record[3]) + int(record[4])
                    depth_control = int(record_control[3]) + int(record_control[4])
                    ratio = float(record[5])
                    ratio_control = float(record_control[5])
                    
                    if depth < self.min_depth or depth_control < self.min_depth or ratio_control == 0:
                        continue
                    
                    # 计算相对位置索引
                    if strand == '+':
                        index = pos - center + self.span // 2
                    else:
                        index = center - pos + self.span // 2
                    
                    if 0 <= index < self.span:
                        self.conversion_list[index].append(ratio / ratio_control)
                except (ValueError, IndexError):
                    continue
        
        print(f"\nProcessed {count} motifs", file=sys.stderr)
        
        # 计算中位数或均值
        if method == 'median':
            self.conversion = np.array([
                np.median(i) if len(i) > 0 else np.nan 
                for i in self.conversion_list
            ])
        else:
            self.conversion = np.array([
                np.mean(i) if len(i) > 0 else np.nan 
                for i in self.conversion_list
            ])
        
        return self.conversion
    
    def plot(self, normalize=False, normalize_flank=False, ylim=None, label=None):
        """
        绘制转化率图
        
        参数：
            normalize: 是否归一化
            normalize_flank: 是否用flank区域归一化
            ylim: Y轴范围
            label: 图例标签
        """
        if self.conversion is None:
            print("Error: No conversion profile calculated", file=sys.stderr)
            return
        
        x = np.arange(-self.span // 2, self.span // 2)
        y = self.conversion.copy()
        
        if normalize:
            if normalize_flank:
                # 使用flank区域归一化（排除中心±20bp）
                flank_indices = (
                    list(range(0, self.span // 2 - 20)) +
                    list(range(self.span // 2 + 20, self.span))
                )
                if len(flank_indices) > 0:
                    factor = np.nanmean(y[flank_indices])
                else:
                    factor = np.nanmean(y)
            else:
                factor = np.nanmean(y)
            
            if not np.isnan(factor) and factor != 0:
                y = y / factor
        
        plt.plot(x, y, label=label if label else 'conversion')
        
        if ylim:
            plt.ylim(ylim)
        
        plt.xlabel('Position relative to motif center (bp)')
        plt.ylabel('Conversion rate' + (' (normalized)' if normalize else ''))
        plt.grid(True, alpha=0.3)
        if label:
            plt.legend()
    
    def close(self):
        """关闭tabix文件"""
        self.gz.close()
        self.gz_control.close()


def main():
    # 检查是否在Snakemake环境中运行
    try:
        # snakemake is available as a global object in Snakemake scripts
        binding_sites = snakemake.input.binding_sites
        sample_files = snakemake.input.samples  # 这是一个列表
        pc_control = snakemake.input.pc_control
        output = snakemake.output[0]
        motif_pick_ratio = snakemake.params.motif_pick_ratio
        min_depth = snakemake.params.min_depth
        span = snakemake.params.span
        normalize_flank = snakemake.params.normalize_flank
        plot_format = snakemake.params.plot_format
        plot_dpi = snakemake.params.plot_dpi
        concentrations = [float(c) for c in snakemake.params.concentrations.split(',')]
        sample_ids = snakemake.params.sample_ids.split(',')
        
        # 确保sample_files是列表
        if not isinstance(sample_files, list):
            sample_files = [sample_files]
    except (NameError, AttributeError):
        # 如果不在Snakemake环境中，使用argparse
        parser = argparse.ArgumentParser(
            description="Plot motif-centered conversion rates at different concentrations"
        )
        parser.add_argument(
            "--binding_sites",
            required=True,
            help="Binding sites BED file (from open region intersection)"
        )
        parser.add_argument(
            "--samples",
            required=True,
            help="Comma-separated list of sample ratio files (tabix compressed)"
        )
        parser.add_argument(
            "--sample_ids",
            required=True,
            help="Comma-separated list of sample IDs"
        )
        parser.add_argument(
            "--concentrations",
            required=True,
            help="Comma-separated list of concentrations (M)"
        )
        parser.add_argument(
            "--pc_control",
            required=True,
            help="PC control ratio file (tabix compressed)"
        )
        parser.add_argument(
            "--output",
            required=True,
            help="Output plot file"
        )
        parser.add_argument(
            "--motif_pick_ratio",
            type=float,
            default=0.01,
            help="Top motif selection ratio"
        )
        parser.add_argument(
            "--min_depth",
            type=int,
            default=5,
            help="Minimum sequencing depth"
        )
        parser.add_argument(
            "--span",
            type=int,
            default=100,
            help="Window span around motif center (bp)"
        )
        parser.add_argument(
            "--normalize_flank",
            action='store_true',
            help="Normalize using flank regions"
        )
        parser.add_argument(
            "--plot_format",
            default="pdf",
            help="Plot format (pdf/png/svg)"
        )
        parser.add_argument(
            "--plot_dpi",
            type=int,
            default=300,
            help="Plot DPI"
        )
        
        args = parser.parse_args()
        binding_sites = args.binding_sites
        sample_files = args.samples.split(',')
        sample_ids = args.sample_ids.split(',')
        concentrations = [float(c) for c in args.concentrations.split(',')]
        pc_control = args.pc_control
        output = args.output
        motif_pick_ratio = args.motif_pick_ratio
        min_depth = args.min_depth
        span = args.span
        normalize_flank = args.normalize_flank
        plot_format = args.plot_format
        plot_dpi = args.plot_dpi
    
    if len(sample_files) != len(sample_ids) or len(sample_files) != len(concentrations):
        print("Error: Number of samples, sample_ids, and concentrations must match", file=sys.stderr)
        sys.exit(1)
    
    # 创建图表
    fig, ax = plt.subplots(figsize=(10, 7))
    plt.sca(ax)
    
    # 为每个浓度绘制转化率曲线
    for sample_file, sample_id, concentration in zip(sample_files, sample_ids, concentrations):
        print(f"\nProcessing sample: {sample_id} (concentration: {concentration} M)", file=sys.stderr)
        
        plotter = MotifCenteredPlotSiteWise(
            sample_file,
            pc_control,
            binding_sites,
            motif_pick_ratio=motif_pick_ratio,
            span=span,
            min_depth=min_depth
        )
        
        try:
            # 选择motifs
            plotter.motifs_list(pick_top_or_bottom='top', motif_pick_ratio=motif_pick_ratio)
            
            # 计算转化率profile
            plotter.conversion_profile(method='median')
            
            # 绘制
            label = f"{sample_id} ({concentration:.2e} M)"
            plotter.plot(
                normalize=True,
                normalize_flank=normalize_flank,
                label=label
            )
        finally:
            plotter.close()
    
    # 设置图表属性
    plt.title(f"Motif-centered conversion rates\n(TF: {Path(binding_sites).parent.parent.parent.name})")
    plt.legend(loc='best')
    plt.tight_layout()
    
    # 保存图表
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output, format=plot_format, dpi=plot_dpi, bbox_inches='tight')
    plt.close()
    
    print(f"\nPlot saved to: {output}", file=sys.stderr)


if __name__ == "__main__":
    main()


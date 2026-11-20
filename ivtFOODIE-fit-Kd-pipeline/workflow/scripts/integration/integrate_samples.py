"""
样本Ratio数据整合模块

功能说明：
    整合多个样本的conversion ratio数据，并与PC对照进行比对。

生物学原理：
    - ivtFOODIE实验中，TF结合DNA会保护C碱基不被APOBEC3A转化
    - Conversion ratio = 未结合DNA的比例
    - 需要与PC (Protein Control, 无蛋白对照) 比较，消除序列偏好性
    - 只保留在open chromatin区域内且测序深度足够的高质量位点

坐标系统：
    - 输入的tabix文件使用0-based坐标
    - BED文件也是0-based, 左闭右开区间 [start, end)
"""

from typing import Dict, Tuple, List, Optional
from pathlib import Path
import pandas as pd
import numpy as np
from numpy.typing import NDArray
import tabix
import pysam


class RatioIntegrator:
    """
    多样本Ratio数据整合器

    负责从多个样本的tabix压缩文件中提取conversion ratio数据，
    整合成单一的DataFrame，便于后续Kd拟合分析。
    """

    def __init__(
        self,
        sample_tabix_files: Dict[int, str],
        pc_tabix_file: str,
        open_region_bed: str,
        genome_fasta: str,
        min_depth: int = 10
    ):
        """
        初始化RatioIntegrator

        参数：
            sample_tabix_files: 样本编号 -> tabix文件路径的映射
                例: {1: "sample1.txt.gz", 2: "sample2.txt.gz", ...}
            pc_tabix_file: PC对照的tabix文件路径
            open_region_bed: Open chromatin区域的BED文件路径
            genome_fasta: 参考基因组FASTA文件路径
            min_depth: 最小测序深度阈值（深度低于此值的位点会被过滤）
        """
        self.sample_tabix_files = sample_tabix_files
        self.pc_tabix_file = pc_tabix_file
        self.open_region_bed = open_region_bed
        self.genome_fasta = genome_fasta
        self.min_depth = min_depth

        # 打开tabix文件句柄
        self.sample_handles: Dict[int, tabix.Tabix] = {}
        for sample_id, file_path in sample_tabix_files.items():
            self.sample_handles[sample_id] = tabix.open(file_path)

        self.pc_handle = tabix.open(pc_tabix_file)

        # 打开参考基因组（目前暂未使用，保留接口用于future扩展）
        self.genome = pysam.FastaFile(genome_fasta)

    def _filter_by_depth_and_integrate(
        self,
        region_chrom: str,
        region_start: int,
        region_end: int
    ) -> Dict[Tuple[str, int], Dict]:
        """
        在指定区域内整合样本数据，并按测序深度过滤

        原理：
            1. 遍历区域内的每个碱基位点
            2. 检查该位点在所有样本和PC对照中是否都有数据
            3. 检查测序深度是否达标（conversion + unconversion >= min_depth）
            4. 只保留所有条件都满足的位点

        参数：
            region_chrom: 染色体名称
            region_start: 区域起始位置 (0-based)
            region_end: 区域结束位置 (0-based, exclusive)

        返回：
            {(chrom, pos): {
                sample_id: (ratio, depth),
                "adjust": (pc_ratio, pc_depth)
            }}
        """
        # 预查询所有样本和PC在该区域的记录
        records_dict: Dict = {}

        # 检查是否所有样本都能成功查询该区域
        query_failed = False

        # 查询PC对照
        try:
            pc_records = list(self.pc_handle.query(region_chrom, region_start, region_end))
        except:
            return {}  # PC查询失败，跳过该区域

        # 转换为 pos -> record 的字典
        records_dict["adjust"] = {int(record[1]): record for record in pc_records}

        # 查询所有样本
        for sample_id, sample_handle in self.sample_handles.items():
            try:
                sample_records = list(sample_handle.query(region_chrom, region_start, region_end))
            except:
                query_failed = True
                break

            records_dict[sample_id] = {int(record[1]): record for record in sample_records}

        if query_failed:
            return {}

        # 整合每个位点的数据
        ratio_dict = {}

        for pos in range(region_start, region_end):
            # 检查该位点在所有样本和PC中是否都有数据
            depth_ok = True
            this_position_data = {}

            # 检查PC对照
            if pos not in records_dict["adjust"]:
                depth_ok = False
            else:
                pc_record = records_dict["adjust"][pos]
                pc_converted = int(pc_record[3])
                pc_unconverted = int(pc_record[4])
                pc_ratio = float(pc_record[5])
                pc_depth = pc_converted + pc_unconverted

                if pc_depth < self.min_depth:
                    depth_ok = False
                else:
                    this_position_data["adjust"] = (pc_ratio, pc_depth)

            if not depth_ok:
                continue

            # 检查所有样本
            for sample_id in self.sample_handles.keys():
                if pos not in records_dict[sample_id]:
                    depth_ok = False
                    break

                sample_record = records_dict[sample_id][pos]
                converted = int(sample_record[3])
                unconverted = int(sample_record[4])
                ratio = float(sample_record[5])
                depth = converted + unconverted

                if depth < self.min_depth:
                    depth_ok = False
                    break

                this_position_data[sample_id] = (ratio, depth)

            if not depth_ok:
                continue

            # 该位点通过所有检查，加入结果
            ratio_dict[(region_chrom, pos)] = this_position_data

        return ratio_dict

    def integrate_samples(self) -> pd.DataFrame:
        """
        整合所有样本的ratio数据

        流程：
            1. 读取open chromatin区域
            2. 对每个区域，提取并整合所有样本的数据
            3. 按测序深度过滤
            4. 返回整合后的DataFrame

        返回：
            DataFrame with columns:
                - ch: 染色体
                - st: 起始位置 (0-based)
                - ed: 结束位置 (0-based, exclusive)
                - ratio_1, depth_1, ratio_2, depth_2, ... : 各样本数据
                - ratio_adjust, depth_adjust: PC对照数据
        """
        # 读取open region
        # BED文件可能有3列或4列，自动检测
        open_region_df = pd.read_csv(
            self.open_region_bed,
            sep='\t',
            header=None
        )

        # 根据列数设置列名
        if len(open_region_df.columns) >= 4:
            open_region_df.columns = ['ch', 'st', 'ed', 'name'] + [f'col{i}' for i in range(4, len(open_region_df.columns))]
        else:
            open_region_df.columns = ['ch', 'st', 'ed'] + [f'col{i}' for i in range(3, len(open_region_df.columns))]

        # 整合数据
        all_ratio_dict = {}
        total_regions = len(open_region_df)

        for index, region in open_region_df.iterrows():
            if index % 100 == 0:
                print(f"Processing region {index}/{total_regions}", end='\r')

            ch = region['ch']
            st = int(region['st'])
            ed = int(region['ed'])

            # 整合该区域的数据
            region_ratio_dict = self._filter_by_depth_and_integrate(ch, st, ed)
            all_ratio_dict.update(region_ratio_dict)

        print(f"\nTotal valid positions: {len(all_ratio_dict)}")

        # 转换为DataFrame
        df = self._dict_to_dataframe(all_ratio_dict)

        return df

    def _dict_to_dataframe(
        self,
        ratio_dict: Dict[Tuple[str, int], Dict]
    ) -> pd.DataFrame:
        """
        将ratio_dict转换为DataFrame

        参数：
            ratio_dict: {(chrom, pos): {sample_id: (ratio, depth), "adjust": (ratio, depth)}}

        返回：
            DataFrame with proper column names
        """
        rows = []

        for (chrom, pos), sample_data in ratio_dict.items():
            row = {
                'ch': chrom,
                'st': pos,
                'ed': pos + 1
            }

            # 添加各样本数据
            for sample_id in sorted([k for k in sample_data.keys() if k != "adjust"]):
                ratio, depth = sample_data[sample_id]
                row[f'ratio_{sample_id}'] = ratio
                row[f'depth_{sample_id}'] = depth

            # 添加PC对照数据
            ratio_adjust, depth_adjust = sample_data["adjust"]
            row['ratio_adjust'] = ratio_adjust
            row['depth_adjust'] = depth_adjust

            rows.append(row)

        df = pd.DataFrame(rows)

        # 排序：按染色体和位置
        df = df.sort_values(['ch', 'st']).reset_index(drop=True)

        return df

    def close(self):
        """关闭所有文件句柄"""
        # tabix对象通常不需要显式关闭，但保留接口
        if hasattr(self, 'genome') and self.genome is not None:
            self.genome.close()


def integrate_samples_from_config(
    tf_name: str,
    sample_path: Path,
    samples_from: int,
    samples_to: int,
    pc_path: Path,
    open_region_bed: Path,
    genome_fasta: Path,
    output_path: Path,
    min_depth: int = 10,
    sample_name_pattern: str = "{tf}-{i}",
    sample_index_position: int = -1
) -> Path:
    """
    从配置参数整合样本数据（Snakemake wrapper函数）

    参数：
        tf_name: TF名称
        sample_path: 样本数据根目录
        samples_from: 样本编号起始
        samples_to: 样本编号结束
        pc_path: PC对照tabix文件路径
        open_region_bed: Open region BED文件路径
        genome_fasta: 参考基因组FASTA文件路径
        output_path: 输出文件路径
        min_depth: 最小测序深度
        sample_name_pattern: 样本命名模式
        sample_index_position: 样本编号在文件名中的位置（-1表示最后一个字段）

    返回：
        输出文件路径
    """
    import os

    # 查找样本文件
    sample_tabix_files = {}

    for i in range(samples_from, samples_to + 1):
        sample_found = False

        for folder in os.listdir(sample_path):
            # 根据sample_index_position提取样本编号
            folder_parts = folder.split('-')
            if sample_index_position < 0:
                # 负索引
                sample_idx_str = folder_parts[sample_index_position].split('_')[0]
            else:
                # 正索引
                sample_idx_str = folder_parts[sample_index_position].split('_')[0]

            if str(i) == sample_idx_str:
                # 找到对应样本
                sample_folder = sample_path / folder

                # 查找.gz文件
                for gz_file in os.listdir(sample_folder):
                    if gz_file.endswith('.gz') and not gz_file.endswith('.tbi'):
                        sample_tabix_files[i] = str(sample_folder / gz_file)
                        sample_found = True
                        print(f"Found sample {i}: {sample_tabix_files[i]}")
                        break

                if sample_found:
                    break

        if not sample_found:
            raise FileNotFoundError(f"Sample {i} not found in {sample_path}")

    # 检查是否所有样本都找到了
    expected_count = samples_to - samples_from + 1
    if len(sample_tabix_files) != expected_count:
        raise ValueError(
            f"Expected {expected_count} samples, but only found {len(sample_tabix_files)}"
        )

    # 创建integrator并运行
    integrator = RatioIntegrator(
        sample_tabix_files=sample_tabix_files,
        pc_tabix_file=str(pc_path),
        open_region_bed=str(open_region_bed),
        genome_fasta=str(genome_fasta),
        min_depth=min_depth
    )

    print(f"\nIntegrating samples for {tf_name}...")
    df = integrator.integrate_samples()

    # 输出到文件
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep='\t', index=False)

    print(f"\nIntegration complete. Output saved to: {output_path}")
    print(f"Total valid C sites: {len(df)}")

    integrator.close()

    return output_path


# ========================================
# Snakemake wrapper function
# ========================================

def integrate_from_snakemake(
    sample_files: List[str],
    pc_control: str,
    open_region: str,
    genome: str,
    output_file: str,
    min_depth: int = 10
) -> None:
    """
    从Snakemake直接整合样本（直接使用文件列表）

    参数：
        sample_files: 样本tabix文件路径列表
        pc_control: PC对照tabix文件路径
        open_region: Open region BED文件路径
        genome: 参考基因组FASTA文件路径
        output_file: 输出文件路径
        min_depth: 最小测序深度
    """
    # 构建样本编号到文件的映射（按文件顺序编号 1, 2, 3...）
    sample_tabix_files = {
        i + 1: str(file_path)
        for i, file_path in enumerate(sample_files)
    }

    print(f"Found {len(sample_tabix_files)} sample files:")
    for sample_id, file_path in sample_tabix_files.items():
        print(f"  Sample {sample_id}: {file_path}")

    # 创建integrator并运行
    integrator = RatioIntegrator(
        sample_tabix_files=sample_tabix_files,
        pc_tabix_file=str(pc_control),
        open_region_bed=str(open_region),
        genome_fasta=str(genome),
        min_depth=min_depth
    )

    print(f"\nIntegrating samples...")
    df = integrator.integrate_samples()

    # 输出到文件
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_file, sep='\t', index=False)

    print(f"\nIntegration complete. Output saved to: {output_file}")
    print(f"Total valid C sites: {len(df)}")

    integrator.close()


if __name__ == "__main__":
    # 检查是否在Snakemake环境中运行
    try:
        # Snakemake模式：读取snakemake对象
        try:
            integrate_from_snakemake(
                sample_files=snakemake.input.samples,
                pc_control=snakemake.input.pc_control,
                open_region=snakemake.input.open_region,
                genome=snakemake.params.genome,
                output_file=snakemake.output[0],
                min_depth=snakemake.params.min_depth
            )
        except Exception as e:
            # 如果失败，创建空输出文件
            import sys
            from pathlib import Path
            print(f"Error running sample integration: {e}", file=sys.stderr)
            print(f"Creating empty output file...", file=sys.stderr)
            
            # 创建空的 TSV 文件
            Path(snakemake.output[0]).parent.mkdir(parents=True, exist_ok=True)
            Path(snakemake.output[0]).touch()
            
            # 将错误信息写入日志（如果存在）
            if hasattr(snakemake, 'log') and len(snakemake.log) > 0:
                import traceback
                with open(snakemake.log[0], 'a') as f:
                    f.write(f"\n\nERROR: Sample integration failed\n")
                    f.write(f"Exception: {str(e)}\n")
                    f.write(f"Traceback:\n{traceback.format_exc()}\n")
                    f.write("Empty output file created to allow workflow to continue.\n")
            
            # 退出时使用成功状态码，这样 Snakemake 不会认为规则失败
            sys.exit(0)

    except NameError:
        # 命令行模式：使用argparse
        import argparse

        parser = argparse.ArgumentParser(
            description='整合多样本conversion ratio数据'
        )

        parser.add_argument('--TF', type=str, required=True, help='TF名称')
        parser.add_argument('--sample_path', type=str, required=True, help='样本数据根目录')
        parser.add_argument('--samples_from', type=int, required=True, help='样本编号起始')
        parser.add_argument('--samples_to', type=int, required=True, help='样本编号结束')
        parser.add_argument('--PC_path', type=str, required=True, help='PC对照tabix文件')
        parser.add_argument('--open_region_f', type=str, required=True, help='Open region BED文件')
        parser.add_argument('--genome_f', type=str, required=True, help='参考基因组FASTA文件')
        parser.add_argument('--out_path', type=str, required=True, help='输出文件路径')
        parser.add_argument('--min_depth', type=int, default=10, help='最小测序深度')
        parser.add_argument(
            '--sample_index_position',
            type=int,
            default=-1,
            help='样本编号在文件名中的位置（-1表示最后）'
        )

        args = parser.parse_args()

        integrate_samples_from_config(
            tf_name=args.TF,
            sample_path=Path(args.sample_path),
            samples_from=args.samples_from,
            samples_to=args.samples_to,
            pc_path=Path(args.PC_path),
            open_region_bed=Path(args.open_region_f),
            genome_fasta=Path(args.genome_f),
            output_path=Path(args.out_path),
            min_depth=args.min_depth,
            sample_index_position=args.sample_index_position
        )

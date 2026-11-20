#!/usr/bin/env python3
"""
提取chr10测试数据脚本

功能：
    1. 从config.yaml读取所有TF配置
    2. 提取每个sample和PC的sites.txt文件的chr10数据（人和小鼠都有chr10）
    3. 保存到test_data目录
    4. 创建测试用的sample sheet
    5. 修改config.yaml指向测试数据
"""

import yaml
import pandas as pd
from pathlib import Path
import argparse
import sys


def extract_chr10_from_sites_file(input_file: Path, output_file: Path, target_chr: str = "chr10") -> int:
    """
    从sites.txt文件中提取指定染色体的数据
    
    参数：
        input_file: 输入的sites.txt文件路径
        output_file: 输出的测试数据文件路径
        target_chr: 目标染色体名称（默认chr10，人和小鼠都有）
    
    返回：
        提取的行数
    """
    print(f"  提取 {input_file} -> {output_file}")
    
    # 读取文件，只提取目标染色体的行
    lines_extracted = 0
    
    # 支持多种染色体命名格式
    target_variants = [
        target_chr,
        target_chr.replace('chr', ''),
        target_chr.lower(),
        target_chr.upper(),
        target_chr.replace('chr', '').lower(),
        target_chr.replace('chr', '').upper()
    ]
    target_variants = list(set(target_variants))  # 去重
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # 检查第一列是否为目标染色体
            parts = line.strip().split('\t')
            if len(parts) > 0:
                chrom = parts[0]
                # 检查是否匹配任何变体
                if chrom in target_variants:
                    outfile.write(line)
                    lines_extracted += 1
    
    print(f"    提取了 {lines_extracted} 行数据")
    return lines_extracted


def get_sample_files_from_config(tf_config: dict, config_dir: Path) -> tuple:
    """
    从TF配置中获取所有sample文件路径和PC文件路径
    
    参数：
        tf_config: TF配置字典
        config_dir: config.yaml所在目录（用于解析相对路径）
    
    返回：
        (sample_files, pc_file): sample文件列表和PC文件路径
    """
    sample_files = []
    pc_file = None
    
    # 获取PC文件
    if "pc_control" in tf_config:
        pc_path = Path(tf_config["pc_control"])
        if not pc_path.is_absolute():
            pc_path = config_dir / pc_path
        pc_file = pc_path
    
    # 获取sample文件
    if "sample_sheet" in tf_config:
        # 方式1：使用sample sheet
        sample_sheet_path = Path(tf_config["sample_sheet"])
        if not sample_sheet_path.is_absolute():
            # 如果路径以config/开头，相对于config_dir的父目录
            if str(sample_sheet_path).startswith("config/"):
                sample_sheet_path = config_dir.parent / sample_sheet_path
            else:
                sample_sheet_path = config_dir / sample_sheet_path
        
        if sample_sheet_path.exists():
            sheet = pd.read_csv(sample_sheet_path, sep='\t')
            base_path = Path(tf_config.get("sample_base_path", ""))
            file_pattern = tf_config.get("sample_file_pattern", "{folder_name}_sites.txt")
            
            for _, row in sheet.iterrows():
                folder_name = row["folder_name"]
                file_path = base_path / file_pattern.format(folder_name=folder_name)
                sample_files.append(file_path)
        else:
            print(f"  警告: Sample sheet不存在: {sample_sheet_path}")
    
    elif "sample_files" in tf_config:
        # 方式2：直接列出文件路径
        for file_path_str in tf_config["sample_files"]:
            file_path = Path(file_path_str)
            if not file_path.is_absolute():
                file_path = config_dir / file_path
            sample_files.append(file_path)
    
    return sample_files, pc_file


def create_test_sample_sheet(tf_name: str, original_sheet: pd.DataFrame, 
                            test_data_dir: Path, output_path: Path) -> pd.DataFrame:
    """
    创建测试用的sample sheet，更新folder_name指向测试数据
    
    参数：
        tf_name: TF名称
        original_sheet: 原始sample sheet的DataFrame
        test_data_dir: 测试数据目录
        output_path: 输出sample sheet路径
    
    返回：
        新的sample sheet DataFrame
    """
    test_sheet = original_sheet.copy()
    
    # 更新folder_name为测试数据路径
    for idx, row in test_sheet.iterrows():
        original_folder = row["folder_name"]
        # 使用sample_id作为测试数据的文件夹名
        test_folder = f"{tf_name}_test_{row['sample_id']}"
        test_sheet.at[idx, "folder_name"] = test_folder
    
    # 保存测试sample sheet
    test_sheet.to_csv(output_path, sep='\t', index=False)
    print(f"  创建测试sample sheet: {output_path}")
    
    return test_sheet


def main():
    parser = argparse.ArgumentParser(
        description='提取chr10测试数据并更新配置'
    )
    parser.add_argument(
        '--config',
        type=str,
        default='config/config.yaml',
        help='配置文件路径（默认: config/config.yaml）'
    )
    parser.add_argument(
        '--target-chr',
        type=str,
        default='chr10',
        help='目标染色体（默认: chr10，人和小鼠都有）'
    )
    parser.add_argument(
        '--test-data-dir',
        type=str,
        default='test_data',
        help='测试数据输出目录（默认: test_data）'
    )
    
    args = parser.parse_args()
    
    # 读取配置文件
    config_path = Path(args.config)
    if not config_path.exists():
        print(f"错误: 配置文件不存在: {config_path}")
        sys.exit(1)
    
    config_dir = config_path.parent
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # 创建测试数据目录
    test_data_dir = Path(args.test_data_dir)
    test_data_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"=" * 60)
    print(f"提取 {args.target_chr} 测试数据（人和小鼠都有此染色体）")
    print(f"=" * 60)
    print(f"配置文件: {config_path}")
    print(f"测试数据目录: {test_data_dir}")
    print()
    
    # 处理每个TF
    tfs = config.get("tfs", {})
    if not tfs:
        print("警告: 配置文件中没有找到TF配置")
        return
    
    for tf_name, tf_config in tfs.items():
        print(f"\n处理 TF: {tf_name}")
        print("-" * 60)
        
        # 获取sample文件和PC文件
        sample_files, pc_file = get_sample_files_from_config(tf_config, config_dir)
        
        # 创建TF的测试数据目录
        tf_test_dir = test_data_dir / tf_name
        tf_test_dir.mkdir(parents=True, exist_ok=True)
        
        # 提取PC文件
        if pc_file and pc_file.exists():
            pc_test_file = tf_test_dir / "PC_control_sites.txt"
            extract_chr10_from_sites_file(pc_file, pc_test_file, args.target_chr)
        elif pc_file:
            print(f"  警告: PC文件不存在: {pc_file}")
        
        # 提取sample文件
        if "sample_sheet" in tf_config:
            # 读取原始sample sheet
            sample_sheet_path = Path(tf_config["sample_sheet"])
            if not sample_sheet_path.is_absolute():
                # 如果路径以config/开头，相对于config_dir的父目录
                if str(sample_sheet_path).startswith("config/"):
                    sample_sheet_path = config_dir.parent / sample_sheet_path
                else:
                    sample_sheet_path = config_dir / sample_sheet_path
            
            if sample_sheet_path.exists():
                original_sheet = pd.read_csv(sample_sheet_path, sep='\t')
                
                # 先创建测试用的sample sheet（确定folder_name）
                # samples目录在config目录下，所以相对于config_dir.parent
                test_sheet_path = config_dir.parent / "config" / "samples" / f"{tf_name}_test_samples.tsv"
                test_sheet_path.parent.mkdir(parents=True, exist_ok=True)
                test_sheet = create_test_sample_sheet(
                    tf_name, original_sheet, tf_test_dir, test_sheet_path
                )
                
                # 提取每个sample的chr10数据（使用测试sample sheet中的folder_name）
                for idx, row in test_sheet.iterrows():
                    if idx < len(sample_files) and sample_files[idx].exists():
                        folder_name = row["folder_name"]
                        test_sample_file = tf_test_dir / f"{folder_name}_sites.txt"
                        extract_chr10_from_sites_file(
                            sample_files[idx], 
                            test_sample_file, 
                            args.target_chr
                        )
                    elif idx < len(sample_files):
                        print(f"  警告: Sample文件不存在: {sample_files[idx]}")
                
                print(f"\n  测试sample sheet已创建: {test_sheet_path}")
                print(f"  共 {len(test_sheet)} 个样本")
            else:
                print(f"  警告: Sample sheet不存在: {sample_sheet_path}")
        
        elif "sample_files" in tf_config:
            # 直接列出文件的方式
            for idx, sample_file in enumerate(sample_files):
                if sample_file.exists():
                    test_sample_file = tf_test_dir / f"sample_{idx+1}_sites.txt"
                    extract_chr10_from_sites_file(
                        sample_file, 
                        test_sample_file, 
                        args.target_chr
                    )
                else:
                    print(f"  警告: Sample文件不存在: {sample_file}")
    
    print(f"\n" + "=" * 60)
    print("提取完成！")
    print("=" * 60)
    print(f"\n测试数据保存在: {test_data_dir}")
    print("\n下一步:")
    print("1. 检查测试数据是否正确提取")
    print("2. 手动更新 config/config.yaml 中的路径指向测试数据")
    print("3. 运行 workflow 进行测试")


if __name__ == "__main__":
    main()


#!/usr/bin/env python3
"""
更新config.yaml以使用测试数据

功能：
    1. 读取config.yaml
    2. 更新TF配置中的路径指向测试数据
    3. 备份原始config.yaml
    4. 保存更新后的config.yaml
"""

import yaml
from pathlib import Path
import argparse
import sys
from datetime import datetime


def update_config_for_test(config_path: Path, test_data_dir: Path, backup: bool = True):
    """
    更新config.yaml以使用测试数据
    
    参数：
        config_path: config.yaml文件路径
        test_data_dir: 测试数据目录
        backup: 是否备份原始配置文件
    """
    # 读取配置文件
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # 备份原始配置
    if backup:
        backup_path = config_path.with_suffix(f'.yaml.backup_{datetime.now().strftime("%Y%m%d_%H%M%S")}')
        with open(backup_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
        print(f"已备份原始配置到: {backup_path}")
    
    # 获取config.yaml的绝对路径（用于计算相对路径）
    config_dir = config_path.parent.absolute()
    test_data_dir_abs = test_data_dir.absolute()
    
    # 计算测试数据目录相对于config.yaml的路径
    try:
        test_data_dir_rel = Path(test_data_dir_abs).relative_to(config_dir.parent)
    except ValueError:
        # 如果无法计算相对路径，使用绝对路径
        test_data_dir_rel = test_data_dir_abs
    
    # 更新每个TF的配置
    tfs = config.get("tfs", {})
    if not tfs:
        print("警告: 配置文件中没有找到TF配置")
        return
    
    for tf_name, tf_config in tfs.items():
        print(f"\n更新 TF: {tf_name}")
        
        # 更新PC文件路径
        if "pc_control" in tf_config:
            original_pc = Path(tf_config["pc_control"])
            test_pc = test_data_dir_abs / tf_name / "PC_control_sites.txt"
            if test_pc.exists():
                # 使用相对路径
                tf_config["pc_control"] = str(test_data_dir_rel / tf_name / "PC_control_sites.txt")
                print(f"  PC文件: {original_pc} -> {tf_config['pc_control']}")
            else:
                print(f"  警告: 测试PC文件不存在: {test_pc}")
        
        # 更新sample sheet路径
        if "sample_sheet" in tf_config:
            original_sheet = Path(tf_config["sample_sheet"])
            test_sheet = config_dir / "samples" / f"{tf_name}_test_samples.tsv"
            if test_sheet.exists():
                tf_config["sample_sheet"] = f"config/samples/{tf_name}_test_samples.tsv"
                print(f"  Sample sheet: {original_sheet} -> {tf_config['sample_sheet']}")
            else:
                print(f"  警告: 测试sample sheet不存在: {test_sheet}")
        
        # 更新sample_base_path
        if "sample_base_path" in tf_config:
            original_base = Path(tf_config["sample_base_path"])
            test_base = test_data_dir_abs / tf_name
            if test_base.exists():
                # 使用相对路径
                tf_config["sample_base_path"] = str(test_data_dir_rel / tf_name)
                print(f"  Sample base path: {original_base} -> {tf_config['sample_base_path']}")
        
        # 更新sample_file_pattern（测试数据使用folder_name命名）
        if "sample_sheet" in tf_config:
            # 测试数据的命名模式：{folder_name}_sites.txt
            # folder_name在测试sample sheet中已经更新为 {tf_name}_test_{sample_id}
            # 所以文件是 {tf_name}_test_{sample_id}_sites.txt
            # pattern应该是 {folder_name}_sites.txt
            tf_config["sample_file_pattern"] = "{folder_name}_sites.txt"
            print(f"  Sample file pattern: {tf_config['sample_file_pattern']}")
    
    # 保存更新后的配置
    with open(config_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
    
    print(f"\n配置已更新: {config_path}")


def main():
    parser = argparse.ArgumentParser(
        description='更新config.yaml以使用测试数据'
    )
    parser.add_argument(
        '--config',
        type=str,
        default='config/config.yaml',
        help='配置文件路径（默认: config/config.yaml）'
    )
    parser.add_argument(
        '--test-data-dir',
        type=str,
        default='test_data',
        help='测试数据目录（默认: test_data）'
    )
    parser.add_argument(
        '--no-backup',
        action='store_true',
        help='不备份原始配置文件'
    )
    
    args = parser.parse_args()
    
    config_path = Path(args.config)
    if not config_path.exists():
        print(f"错误: 配置文件不存在: {config_path}")
        sys.exit(1)
    
    test_data_dir = Path(args.test_data_dir)
    if not test_data_dir.exists():
        print(f"错误: 测试数据目录不存在: {test_data_dir}")
        print("请先运行 extract_chr20_test_data.py (或 extract_chr10_test_data.py) 提取测试数据")
        sys.exit(1)
    
    print("=" * 60)
    print("更新config.yaml以使用测试数据")
    print("=" * 60)
    print(f"配置文件: {config_path}")
    print(f"测试数据目录: {test_data_dir}")
    print()
    
    update_config_for_test(config_path, test_data_dir, backup=not args.no_backup)
    
    print("\n" + "=" * 60)
    print("更新完成！")
    print("=" * 60)


if __name__ == "__main__":
    main()


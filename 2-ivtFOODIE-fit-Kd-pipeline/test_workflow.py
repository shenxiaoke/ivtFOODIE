#!/usr/bin/env python3
"""
测试workflow配置和辅助函数

用法：
    python test_workflow.py
"""

import sys
from pathlib import Path

# 添加workflow脚本路径
workflow_dir = Path(__file__).parent / "workflow"
sys.path.insert(0, str(workflow_dir / "scripts" / "utils"))

import yaml
import pandas as pd


def test_config_loading():
    """测试配置文件加载"""
    print("="*60)
    print("测试1: 配置文件加载")
    print("="*60)

    config_path = Path("config/config.yaml")
    with open(config_path) as f:
        config = yaml.safe_load(f)

    print(f"✓ 配置文件加载成功")
    print(f"  - TFs配置数量: {len(config['tfs'])}")
    print(f"  - TF名称: {list(config['tfs'].keys())}")
    print()


def test_sample_sheet():
    """测试sample sheet加载"""
    print("="*60)
    print("测试2: Sample Sheet加载")
    print("="*60)

    config_path = Path("config/config.yaml")
    with open(config_path) as f:
        config = yaml.safe_load(f)

    for tf_name, tf_config in config['tfs'].items():
        print(f"\nTF: {tf_name}")

        if "sample_sheet" in tf_config:
            sheet_path = tf_config["sample_sheet"]
            print(f"  使用Sample Sheet: {sheet_path}")

            sheet = pd.read_csv(sheet_path, sep='\t')
            print(f"  样本数量: {len(sheet)}")
            print(f"  列: {list(sheet.columns)}")
            print(f"\n  样本列表:")
            for idx, row in sheet.iterrows():
                print(f"    [{idx}] {row['sample_id']}: {row['concentration']} ({row['concentration_M']:.2e} M)")

            # 测试文件路径生成
            base_path = Path(tf_config.get("sample_base_path", ""))
            file_pattern = tf_config.get("sample_file_pattern", "{folder_name}_sites.txt")

            print(f"\n  文件路径测试:")
            print(f"    Base path: {base_path}")
            print(f"    Pattern: {file_pattern}")

            for idx, row in sheet.iterrows():
                folder_name = row["folder_name"]
                file_path = base_path / file_pattern.format(folder_name=folder_name)
                exists = "✓" if file_path.exists() else "✗"
                print(f"    {exists} {file_path}")
                if idx >= 2:  # 只显示前3个
                    remaining = len(sheet) - 3
                    if remaining > 0:
                        print(f"    ... 还有 {remaining} 个文件")
                    break

    print()


def test_sample_info_function():
    """测试get_sample_info函数"""
    print("="*60)
    print("测试3: get_sample_info函数")
    print("="*60)

    # 加载配置
    config_path = Path("config/config.yaml")
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # 模拟Snakefile中的get_sample_info函数
    def load_sample_sheet(tf_name):
        tf_config = config["tfs"][tf_name]
        sample_sheet_path = tf_config.get("sample_sheet")
        if not sample_sheet_path:
            return None
        sheet = pd.read_csv(sample_sheet_path, sep='\t')
        return sheet

    def get_sample_info(tf_name):
        tf_config = config["tfs"][tf_name]

        if "sample_sheet" in tf_config:
            sheet = load_sample_sheet(tf_name)
            base_path = Path(tf_config.get("sample_base_path", ""))
            file_pattern = tf_config.get("sample_file_pattern", "{folder_name}_sites.txt")

            sample_files = []
            for _, row in sheet.iterrows():
                folder_name = row["folder_name"]
                file_path = base_path / file_pattern.format(folder_name=folder_name)
                sample_files.append(str(file_path))

            return {
                'sample_ids': sheet["sample_id"].tolist(),
                'sample_files': sample_files,
                'concentrations': sheet["concentration_M"].tolist(),
                'fitting_indices': tf_config.get("fitting_indices", list(range(len(sheet))))
            }
        else:
            raise ValueError(f"TF '{tf_name}': Must specify 'sample_sheet'")

    # 测试
    for tf_name in config['tfs'].keys():
        print(f"\nTF: {tf_name}")
        info = get_sample_info(tf_name)

        print(f"  Sample IDs: {info['sample_ids']}")
        print(f"  浓度数量: {len(info['concentrations'])}")
        print(f"  浓度范围: {min(info['concentrations']):.2e} - {max(info['concentrations']):.2e} M")
        print(f"  Fitting indices: {info['fitting_indices']}")
        print(f"  文件数量: {len(info['sample_files'])}")

    print()


def test_config_validation():
    """测试配置验证"""
    print("="*60)
    print("测试4: 配置验证")
    print("="*60)

    from validate_config import validate_config

    config_path = Path("config/config.yaml")
    with open(config_path) as f:
        config = yaml.safe_load(f)

    valid = validate_config(config)

    if valid:
        print("\n✅ 所有配置验证通过！")
    else:
        print("\n❌ 配置验证失败！")

    print()


if __name__ == "__main__":
    try:
        test_config_loading()
        test_sample_sheet()
        test_sample_info_function()
        test_config_validation()

        print("="*60)
        print("✅ 所有测试通过！")
        print("="*60)
        print("\n现在可以运行: snakemake -n")

    except Exception as e:
        print(f"\n❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

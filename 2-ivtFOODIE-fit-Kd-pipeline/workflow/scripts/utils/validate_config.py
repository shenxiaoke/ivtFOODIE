"""
配置文件校验模块

功能：
    在Snakemake workflow启动前校验所有关键文件路径的存在性
    避免在运行到后期阶段才发现路径错误
"""

from pathlib import Path
from typing import Dict, List, Tuple
import sys


class ConfigValidator:
    """配置文件校验器"""

    def __init__(self, config: Dict, base_dir: Path = None):
        """
        初始化配置校验器

        参数：
            config: Snakemake配置字典
            base_dir: 基准目录，用于解析相对路径（默认：当前工作目录）
        """
        self.config = config
        self.base_dir = Path(base_dir) if base_dir else Path.cwd()
        self.errors: List[str] = []
        self.warnings: List[str] = []

    def _check_file_exists(self, file_path: str, description: str) -> bool:
        """
        检查文件是否存在

        参数：
            file_path: 文件路径（可以是相对路径或绝对路径）
            description: 文件描述（用于错误消息）

        返回：
            文件是否存在
        """
        # 解析路径：如果是相对路径，相对于base_dir
        path = Path(file_path)
        if not path.is_absolute():
            path = self.base_dir / path
        
        if not path.exists():
            self.errors.append(f"{description} not found: {file_path}")
            return False
        return True

    def _check_dir_exists(self, dir_path: str, description: str) -> bool:
        """
        检查目录是否存在

        参数：
            dir_path: 目录路径
            description: 目录描述

        返回：
            目录是否存在
        """
        if not Path(dir_path).exists():
            self.errors.append(f"{description} not found: {dir_path}")
            return False
        if not Path(dir_path).is_dir():
            self.errors.append(f"{description} is not a directory: {dir_path}")
            return False
        return True

    def _check_tabix_index_if_gz(self, file_path: str, description: str) -> bool:
        """检查已压缩ratio文件是否带有对应的tabix索引"""
        if not str(file_path).endswith(".gz"):
            return True
        return self._check_file_exists(f"{file_path}.tbi", f"{description} tabix index")

    def validate_genomes(self) -> bool:
        """
        校验参考基因组文件

        返回：
            所有基因组文件是否都存在
        """
        print("\n[1/4] Validating genome files...")

        all_valid = True
        genomes = self.config.get("genomes", {})

        if not genomes:
            self.warnings.append("No genome files specified in config")
            return True

        for species, genome_path in genomes.items():
            if not self._check_file_exists(genome_path, f"Genome file for {species}"):
                all_valid = False
            else:
                # 检查索引文件(.fai)
                fai_path = f"{genome_path}.fai"
                if not Path(fai_path).exists():
                    self.warnings.append(
                        f"Genome index (.fai) not found for {species}: {fai_path}\n"
                        f"  You may need to run: samtools faidx {genome_path}"
                    )

        return all_valid

    def validate_open_regions(self) -> bool:
        """
        校验open chromatin区域文件

        返回：
            所有open region文件是否都存在
        """
        print("[2/4] Validating open region files...")

        all_valid = True
        open_regions = self.config.get("open_regions", {})

        if not open_regions:
            self.warnings.append("No open region files specified in config")
            return True

        # 收集TF配置中使用的cell lines
        used_cell_lines = set()
        for tf_name, tf_config in self.config.get("tfs", {}).items():
            cell_line = tf_config.get("cell_line")
            if cell_line:
                used_cell_lines.add(cell_line)

        # 只检查被使用的cell lines
        for cell_line in used_cell_lines:
            if cell_line not in open_regions:
                self.errors.append(
                    f"Open region file for cell line '{cell_line}' not defined in config"
                )
                all_valid = False
            else:
                region_path = open_regions[cell_line]
                if not self._check_file_exists(
                    region_path,
                    f"Open region file for {cell_line}"
                ):
                    all_valid = False

        return all_valid

    def validate_tf_configs(self) -> bool:
        """
        校验TF特异性配置

        返回：
            所有TF配置是否有效
        """
        print("[3/4] Validating TF configurations...")

        all_valid = True
        tfs = self.config.get("tfs", {})

        if not tfs:
            self.errors.append("No TFs defined in config")
            return False

        for tf_name, tf_config in tfs.items():
            print(f"  Checking TF: {tf_name}")

            # 检查基本必需字段
            basic_required = ["species", "cell_line", "pc_control"]
            for field in basic_required:
                if field not in tf_config:
                    self.errors.append(f"TF '{tf_name}': missing required field '{field}'")
                    all_valid = False

            # 检查PC对照文件（txt.gz输入必须带有预先建立的tabix索引）
            pc_control = tf_config.get("pc_control")
            if pc_control:
                pc_description = f"PC control file for {tf_name}"
                if not self._check_file_exists(pc_control, pc_description):
                    all_valid = False
                elif not self._check_tabix_index_if_gz(pc_control, pc_description):
                    all_valid = False

            # 方式1：使用sample sheet
            if "sample_sheet" in tf_config:
                sample_sheet = tf_config["sample_sheet"]

                # 解析sample sheet路径
                sample_sheet_path = Path(sample_sheet)
                if not sample_sheet_path.is_absolute():
                    sample_sheet_path = self.base_dir / sample_sheet_path

                # 检查sample sheet文件是否存在
                if not self._check_file_exists(sample_sheet, f"Sample sheet for {tf_name}"):
                    all_valid = False
                else:
                    # 验证sample sheet格式
                    try:
                        import pandas as pd
                        sheet = pd.read_csv(sample_sheet_path, sep='\t')

                        # 检查必需列
                        required_columns = ["sample_id", "concentration_M", "folder_name"]
                        for col in required_columns:
                            if col not in sheet.columns:
                                self.errors.append(
                                    f"TF '{tf_name}': Sample sheet missing column '{col}'"
                                )
                                all_valid = False

                        # 检查base_path
                        base_path = tf_config.get("sample_base_path")
                        if not base_path:
                            self.warnings.append(
                                f"TF '{tf_name}': 'sample_base_path' not specified, using relative paths"
                            )
                        elif not Path(base_path).exists():
                            self.errors.append(
                                f"TF '{tf_name}': sample_base_path not found: {base_path}"
                            )
                            all_valid = False
                        else:
                            # 检查每个样本文件是否存在
                            file_pattern = tf_config.get("sample_file_pattern", "{folder_name}_sites.txt")
                            for idx, row in sheet.iterrows():
                                folder_name = row["folder_name"]
                                file_path = Path(base_path) / file_pattern.format(folder_name=folder_name)
                                if not file_path.exists():
                                    self.errors.append(
                                        f"TF '{tf_name}': Sample file not found: {file_path}"
                                    )
                                    all_valid = False
                                elif not self._check_tabix_index_if_gz(
                                    str(file_path),
                                    f"TF '{tf_name}': Sample file"
                                ):
                                    all_valid = False

                        # 检查fitting_indices
                        fitting_indices = tf_config.get("fitting_indices", [])
                        if fitting_indices:
                            max_index = max(fitting_indices)
                            if max_index >= len(sheet):
                                self.errors.append(
                                    f"TF '{tf_name}': fitting_indices contains index {max_index}, "
                                    f"but sample sheet only has {len(sheet)} samples (0-indexed)"
                                )
                                all_valid = False
                        else:
                            self.warnings.append(
                                f"TF '{tf_name}': 'fitting_indices' not specified, will use all samples"
                            )

                    except Exception as e:
                        self.errors.append(
                            f"TF '{tf_name}': Error reading sample sheet: {e}"
                        )
                        all_valid = False

            # 方式2：直接列出sample_files
            elif "sample_files" in tf_config:
                sample_files = tf_config.get("sample_files")

                if not isinstance(sample_files, list):
                    self.errors.append(
                        f"TF '{tf_name}': 'sample_files' must be a list"
                    )
                    all_valid = False
                else:
                    # 检查每个样本文件是否存在
                    for i, sample_file in enumerate(sample_files, 1):
                        sample_description = f"Sample {i} file for {tf_name}"
                        if not self._check_file_exists(
                            sample_file,
                            sample_description
                        ):
                            all_valid = False
                        elif not self._check_tabix_index_if_gz(sample_file, sample_description):
                            all_valid = False

                    # 检查concentrations
                    concentrations = tf_config.get("concentrations", [])
                    if not concentrations:
                        self.errors.append(
                            f"TF '{tf_name}': 'concentrations' is required when using 'sample_files'"
                        )
                        all_valid = False
                    elif len(concentrations) != len(sample_files):
                        self.errors.append(
                            f"TF '{tf_name}': concentration count ({len(concentrations)}) does not match "
                            f"sample_files count ({len(sample_files)})"
                        )
                        all_valid = False

                    # 检查fitting_indices
                    fitting_indices = tf_config.get("fitting_indices", [])
                    if fitting_indices:
                        max_index = max(fitting_indices)
                        if max_index >= len(sample_files):
                            self.errors.append(
                                f"TF '{tf_name}': fitting_indices contains index {max_index}, "
                                f"but only {len(sample_files)} samples available (0-indexed)"
                            )
                            all_valid = False

            else:
                # 既没有sample_sheet也没有sample_files
                self.errors.append(
                    f"TF '{tf_name}': Must specify either 'sample_sheet' or 'sample_files'"
                )
                all_valid = False

            # 检查preliminary motif文件（如果指定）
            preliminary_motif = tf_config.get("preliminary_motif")
            if preliminary_motif and preliminary_motif != "null":
                if not self._check_file_exists(
                    preliminary_motif,
                    f"Preliminary motif file for {tf_name}"
                ):
                    all_valid = False

            # 检查species是否在genomes中定义
            species = tf_config.get("species")
            if species and species not in self.config.get("genomes", {}):
                self.errors.append(
                    f"TF '{tf_name}': species '{species}' not defined in genomes section"
                )
                all_valid = False

            # 检查cell_line是否在open_regions中定义
            cell_line = tf_config.get("cell_line")
            if cell_line and cell_line not in self.config.get("open_regions", {}):
                self.errors.append(
                    f"TF '{tf_name}': cell_line '{cell_line}' not defined in open_regions section"
                )
                all_valid = False

        return all_valid

    def validate_output_dirs(self) -> bool:
        """
        检查输出目录是否可写

        返回：
            输出目录是否可访问
        """
        print("[4/4] Validating output directories...")

        all_valid = True
        output_dirs = ["results", "logs", "benchmarks"]

        for dir_name in output_dirs:
            dir_path = Path(dir_name)

            # 如果目录不存在，尝试创建
            if not dir_path.exists():
                try:
                    dir_path.mkdir(parents=True, exist_ok=True)
                    print(f"  Created directory: {dir_path}")
                except Exception as e:
                    self.errors.append(f"Cannot create output directory {dir_path}: {e}")
                    all_valid = False

            # 检查是否可写
            elif not dir_path.is_dir():
                self.errors.append(f"Output path exists but is not a directory: {dir_path}")
                all_valid = False

        return all_valid

    def validate_all(self) -> bool:
        """
        执行所有校验

        返回：
            所有校验是否通过
        """
        print("\n" + "="*60)
        print("Validating configuration...")
        print("="*60)

        results = [
            self.validate_genomes(),
            self.validate_open_regions(),
            self.validate_tf_configs(),
            self.validate_output_dirs()
        ]

        all_valid = all(results)

        # 打印结果
        print("\n" + "="*60)
        if self.warnings:
            print(f"\n⚠️  {len(self.warnings)} Warning(s):")
            for warning in self.warnings:
                print(f"  - {warning}")

        if self.errors:
            print(f"\n❌ {len(self.errors)} Error(s):")
            for error in self.errors:
                print(f"  - {error}")
            print("\n" + "="*60)
            print("Configuration validation FAILED")
            print("Please fix the errors above before running the workflow")
            print("="*60 + "\n")
        else:
            print("\n✅ Configuration validation PASSED")
            if self.warnings:
                print("⚠️  Please review warnings above")
            print("="*60 + "\n")

        return all_valid


def validate_config(config: Dict, base_dir: Path = None) -> bool:
    """
    校验配置文件（便捷函数）

    参数：
        config: Snakemake配置字典
        base_dir: 基准目录，用于解析相对路径（默认：当前工作目录）

    返回：
        校验是否通过
    """
    validator = ConfigValidator(config, base_dir)
    return validator.validate_all()


if __name__ == "__main__":
    # 命令行模式：从YAML文件读取配置并校验
    import yaml
    import argparse

    parser = argparse.ArgumentParser(description="Validate ivtFOODIE workflow configuration")
    parser.add_argument(
        "--config",
        type=str,
        default="config/config.yaml",
        help="Path to config.yaml file"
    )

    args = parser.parse_args()

    # 读取配置
    config_path = Path(args.config)
    if not config_path.exists():
        print(f"Error: Config file not found: {config_path}")
        sys.exit(1)

    with open(config_path) as f:
        config = yaml.safe_load(f)

    # 校验
    valid = validate_config(config)

    # 返回状态码
    sys.exit(0 if valid else 1)

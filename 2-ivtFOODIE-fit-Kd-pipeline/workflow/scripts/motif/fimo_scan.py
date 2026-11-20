"""
FIMO motif scanning script

功能：
    使用FIMO在基因组上扫描motif结合位点，
    并将结果转换为BED格式
"""

import argparse
import subprocess
import sys
from pathlib import Path
import pandas as pd
import re


def extract_motifs_from_meme(meme_file: str):
    """
    从MEME文件中提取所有motif名称
    
    参数：
        meme_file: MEME格式文件路径
    
    返回：
        motif名称列表
    """
    motifs = []
    with open(meme_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("MOTIF"):
                # 格式: "MOTIF motif_name" 或 "MOTIF motif_name description"
                parts = line.split()
                if len(parts) > 1:
                    motif_name = parts[1]
                    motifs.append(motif_name)
    return motifs


def check_motif_exists(meme_file: str, motif_name: str) -> bool:
    """
    检查motif是否存在于MEME文件中
    
    参数：
        meme_file: 输入MEME文件
        motif_name: 要检查的motif名称
    
    返回：
        如果motif存在返回True，否则返回False
    """
    try:
        with open(meme_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith("MOTIF"):
                    parts = line.split()
                    if len(parts) > 1 and parts[1] == motif_name:
                        return True
    except FileNotFoundError:
        print(f"Error: MEME file not found: {meme_file}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"Error reading MEME file {meme_file}: {e}", file=sys.stderr)
        return False
    
    return False


def get_available_motifs(meme_file: str) -> list:
    """
    获取MEME文件中所有可用的motif名称
    
    参数：
        meme_file: 输入MEME文件
    
    返回：
        motif名称列表
    """
    motifs = []
    try:
        with open(meme_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith("MOTIF"):
                    parts = line.split()
                    if len(parts) > 1:
                        motifs.append(parts[1])
    except Exception as e:
        print(f"Error reading MEME file {meme_file}: {e}", file=sys.stderr)
    
    return motifs


def extract_single_motif_from_meme(meme_file: str, motif_name: str, output_meme: str):
    """
    从MEME文件中提取单个motif并创建新的MEME文件
    
    参数：
        meme_file: 输入MEME文件
        motif_name: 要提取的motif名称
        output_meme: 输出MEME文件路径
    """
    # 先验证motif是否存在
    if not check_motif_exists(meme_file, motif_name):
        available_motifs = get_available_motifs(meme_file)
        print(f"Error: Motif '{motif_name}' not found in {meme_file}", file=sys.stderr)
        if available_motifs:
            print(f"Available motifs in {meme_file}: {', '.join(available_motifs)}", file=sys.stderr)
        else:
            print(f"Warning: No motifs found in {meme_file} or file is empty", file=sys.stderr)
        sys.exit(1)
    
    # 先读取整个文件
    with open(meme_file, 'r') as f:
        lines = f.readlines()
    
    # 找到头部结束位置（第一个MOTIF之前）
    header_end = 0
    for i, line in enumerate(lines):
        if line.startswith("MOTIF"):
            header_end = i
            break
    
    # 找到目标motif的位置
    motif_start = None
    motif_end = None
    
    for i, line in enumerate(lines):
        if line.startswith("MOTIF"):
            parts = line.split()
            if len(parts) > 1 and parts[1] == motif_name:
                motif_start = i
                # 找到下一个MOTIF或文件结束
                for j in range(i + 1, len(lines)):
                    if lines[j].startswith("MOTIF"):
                        motif_end = j
                        break
                if motif_end is None:
                    motif_end = len(lines)
                break
    
    if motif_start is None:
        # 这不应该发生，因为我们已经检查过了，但为了安全起见还是保留这个检查
        available_motifs = get_available_motifs(meme_file)
        print(f"Error: Motif '{motif_name}' not found in {meme_file}", file=sys.stderr)
        if available_motifs:
            print(f"Available motifs in {meme_file}: {', '.join(available_motifs)}", file=sys.stderr)
        sys.exit(1)
    
    # 写入输出文件
    with open(output_meme, 'w') as f_out:
        # 写入头部
        f_out.writelines(lines[:header_end])
        # 写入目标motif
        f_out.writelines(lines[motif_start:motif_end])


def run_fimo(
    meme_file: str,
    genome_fasta: str,
    output_dir: str,
    thresh: float = 0.0001,
    max_stored: int = 100000000
):
    """
    运行FIMO进行基因组扫描
    
    参数：
        meme_file: 输入MEME格式motif文件
        genome_fasta: 基因组FASTA文件
        output_dir: FIMO输出目录
        thresh: p-value阈值
        max_stored: 最大存储结果数
    """
    cmd = [
        "fimo",
        "--thresh", str(thresh),
        "--max-stored-scores", str(max_stored),
        "--verbosity", "5",
        "--oc", output_dir,
        meme_file,
        genome_fasta
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running FIMO: {result.stderr}", file=sys.stderr)
        sys.exit(1)
    
    print(f"FIMO completed. Output directory: {output_dir}")


def fimo_to_bed(fimo_dir: str, output_bed: str, motif_name: str = None):
    """
    将FIMO输出转换为BED格式
    
    FIMO输出文件：
        - fimo.tsv: TSV格式的结果文件
    
    BED格式：
        chr, start, end, motif_name, score, strand
    
    参数：
        fimo_dir: FIMO输出目录
        output_bed: 输出BED文件路径
        motif_name: 如果指定，只提取该motif的结果
    """
    fimo_tsv = Path(fimo_dir) / "fimo.tsv"
    
    if not fimo_tsv.exists():
        print(f"Warning: {fimo_tsv} not found. Creating empty BED file.", file=sys.stderr)
        # 创建空的BED文件
        Path(output_bed).parent.mkdir(parents=True, exist_ok=True)
        with open(output_bed, 'w') as f:
            pass
        return
    
    # 读取FIMO TSV文件
    # FIMO TSV格式：
    # motif_id, motif_alt_id, sequence_name, start, stop, strand, score, p-value, q-value, matched_sequence
    try:
        # 先检查文件是否为空或只有注释行
        with open(fimo_tsv, 'r') as f:
            has_data = False
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    has_data = True
                    break
        
        if not has_data:
            print(f"Warning: FIMO output file is empty or contains only comments. Creating empty BED file.", file=sys.stderr)
            Path(output_bed).parent.mkdir(parents=True, exist_ok=True)
            with open(output_bed, 'w') as f:
                pass
            return
        
        # 尝试读取TSV文件
        df = pd.read_csv(fimo_tsv, sep='\t', comment='#')
        
    except (pd.errors.EmptyDataError, pd.errors.ParserError) as e:
        # 如果文件为空或无法解析，创建空BED文件
        print(f"Warning: FIMO TSV file is empty or cannot be parsed: {e}. Creating empty BED file.", file=sys.stderr)
        Path(output_bed).parent.mkdir(parents=True, exist_ok=True)
        with open(output_bed, 'w') as f:
            pass
        return
    except Exception as e:
        print(f"Error reading FIMO TSV: {e}", file=sys.stderr)
        sys.exit(1)
    
    if len(df) == 0 or df.empty:
        print(f"Warning: No results in FIMO output. Creating empty BED file.", file=sys.stderr)
        Path(output_bed).parent.mkdir(parents=True, exist_ok=True)
        with open(output_bed, 'w') as f:
            pass
        return
    
    # 检查必需的列
    required_cols = ['sequence_name', 'start', 'stop', 'strand', 'score']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Error: Missing required columns in FIMO output: {missing_cols}", file=sys.stderr)
        print(f"Available columns: {df.columns.tolist()}", file=sys.stderr)
        sys.exit(1)
    
    # 获取motif名称列（可能是motif_id或motif_alt_id）
    if 'motif_id' in df.columns:
        motif_col = 'motif_id'
    elif 'motif_alt_id' in df.columns:
        motif_col = 'motif_alt_id'
    else:
        # 如果没有motif名称列，使用默认值
        df['motif_name'] = 'motif_1'
        motif_col = 'motif_name'
    
    # 如果指定了motif_name，只提取该motif的结果
    if motif_name is not None:
        df = df[df[motif_col] == motif_name].copy()
        if len(df) == 0:
            print(f"Warning: No results for motif '{motif_name}'. Creating empty BED file.", file=sys.stderr)
            Path(output_bed).parent.mkdir(parents=True, exist_ok=True)
            with open(output_bed, 'w') as f:
                pass
            return
    
    # 构建BED格式DataFrame
    bed_df = pd.DataFrame({
        'chr': df['sequence_name'],
        'start': df['start'] - 1,  # BED是0-based，FIMO是1-based
        'end': df['stop'],
        'motif_name': df[motif_col] if motif_name is None else motif_name,
        'score': df['score'],
        'strand': df['strand']
    })
    
    # 如果有matched_sequence列，可以添加（作为第7列）
    if 'matched_sequence' in df.columns:
        bed_df['sequence'] = df['matched_sequence']
    
    # 按染色体和位置排序
    bed_df = bed_df.sort_values(['chr', 'start', 'end'])
    
    # 写入BED文件
    Path(output_bed).parent.mkdir(parents=True, exist_ok=True)
    bed_df.to_csv(output_bed, sep='\t', header=False, index=False)
    
    print(f"Converted {len(bed_df)} binding sites to BED format: {output_bed}")


def main():
    # 检查是否在Snakemake环境中运行
    try:
        # snakemake is available as a global object in Snakemake scripts
        meme_file = snakemake.input.meme
        genome_fasta = snakemake.params.genome
        output_bed = snakemake.output.bed
        fimo_dir = str(snakemake.output.fimo_dir)
        thresh = snakemake.params.thresh
        max_stored = snakemake.params.max_stored
        # 检查是否有motif_name参数（用于单个motif扫描）
        motif_name = getattr(snakemake.params, 'motif_name', None)
    except (NameError, AttributeError):
        # 如果不在Snakemake环境中，使用argparse
        parser = argparse.ArgumentParser(
            description="Run FIMO motif scanning on genome"
        )
        parser.add_argument(
            "--meme_file",
            required=True,
            help="Input MEME format motif file"
        )
        parser.add_argument(
            "--genome_fasta",
            required=True,
            help="Genome FASTA file"
        )
        parser.add_argument(
            "--output_bed",
            required=True,
            help="Output BED file"
        )
        parser.add_argument(
            "--fimo_dir",
            required=True,
            help="FIMO output directory"
        )
        parser.add_argument(
            "--thresh",
            type=float,
            default=0.0001,
            help="P-value threshold"
        )
        parser.add_argument(
            "--max_stored",
            type=int,
            default=100000000,
            help="Maximum stored scores"
        )
        parser.add_argument(
            "--motif_name",
            type=str,
            required=True,
            help="Motif name to scan (required)"
        )
        
        args = parser.parse_args()
        meme_file = args.meme_file
        genome_fasta = args.genome_fasta
        output_bed = args.output_bed
        fimo_dir = args.fimo_dir
        thresh = args.thresh
        max_stored = args.max_stored
        motif_name = args.motif_name
    
    # 创建输出目录
    Path(fimo_dir).mkdir(parents=True, exist_ok=True)
    Path(output_bed).parent.mkdir(parents=True, exist_ok=True)
    
    # motif_name 是必需的参数（此脚本只处理单个motif）
    if motif_name is None:
        print("Error: motif_name parameter is required. This script only processes single motifs.", file=sys.stderr)
        print("Use fimo_scan_per_motif rule for processing individual motifs.", file=sys.stderr)
        sys.exit(1)
    
    # 验证MEME文件存在
    if not Path(meme_file).exists():
        print(f"Error: MEME file not found: {meme_file}", file=sys.stderr)
        sys.exit(1)
    
    # 在提取motif之前先验证motif是否存在（提前失败，避免浪费资源）
    print(f"Checking if motif '{motif_name}' exists in {meme_file}...")
    if not check_motif_exists(meme_file, motif_name):
        available_motifs = get_available_motifs(meme_file)
        print(f"Error: Motif '{motif_name}' not found in {meme_file}", file=sys.stderr)
        if available_motifs:
            print(f"Available motifs in {meme_file}: {', '.join(available_motifs)}", file=sys.stderr)
        else:
            print(f"Warning: No motifs found in {meme_file} or file is empty", file=sys.stderr)
        print(f"This usually indicates a mismatch between the motif list and the MEME file.", file=sys.stderr)
        print(f"Please check the aggregate_motif_outputs rule to ensure it extracts motifs from the correct MEME file.", file=sys.stderr)
        sys.exit(1)
    
    # 处理单个motif
    temp_meme = str(Path(fimo_dir).parent / f"temp_{motif_name}.meme")
    print(f"Extracting motif '{motif_name}' from MEME file...")
    extract_single_motif_from_meme(meme_file, motif_name, temp_meme)
    meme_file_to_use = temp_meme
    
    # 步骤1: 运行FIMO
    print("Running FIMO...")
    run_fimo(
        meme_file_to_use,
        genome_fasta,
        fimo_dir,
        thresh=thresh,
        max_stored=max_stored
    )
    
    # 步骤2: 转换为BED格式
    print("Converting FIMO output to BED format...")
    fimo_to_bed(fimo_dir, output_bed, motif_name=motif_name)
    
    # 清理临时文件
    if Path(temp_meme).exists():
        Path(temp_meme).unlink()
    
    print("FIMO scanning completed successfully!")


if __name__ == "__main__":
    main()


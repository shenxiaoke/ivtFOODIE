"""
STREME motif enrichment script

功能：
    从footprint BED文件中提取序列，使用STREME进行motif富集，
    并将结果转换为MEME格式
"""

import argparse
import subprocess
import sys
from pathlib import Path
import pandas as pd


def extract_sequences(bed_file: str, genome_fasta: str, output_fasta: str):
    """
    使用bedtools getfasta从BED文件提取序列
    
    参数：
        bed_file: 输入BED文件路径
        genome_fasta: 基因组FASTA文件路径
        output_fasta: 输出FASTA文件路径
    """
    cmd = [
        "bedtools", "getfasta",
        "-fi", genome_fasta,
        "-bed", bed_file,
        "-fo", output_fasta,
        "-name"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running bedtools getfasta: {result.stderr}", file=sys.stderr)
        sys.exit(1)


def run_streme(
    input_fasta: str,
    output_dir: str,
    minw: int = 6,
    maxw: int = 15,
    nmotifs: int = 10,
    thresh: float = 0.05,
    totallength: int = 4000000
):
    """
    运行STREME进行motif富集
    
    参数：
        input_fasta: 输入FASTA文件（正序列）
        output_dir: STREME输出目录
        minw: 最小motif宽度
        maxw: 最大motif宽度
        nmotifs: 发现motif数量
        thresh: p-value阈值
        totallength: 用于计算p-value的总序列长度
    """
    cmd = [
        "streme",
        "--p", input_fasta,
        "--oc", output_dir,
        "--minw", str(minw),
        "--maxw", str(maxw),
        "--nmotifs", str(nmotifs),
        "--thresh", str(thresh),
        "--totallength", str(totallength)
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running STREME: {result.stderr}", file=sys.stderr)
        sys.exit(1)
    
    print(f"STREME completed. Output directory: {output_dir}")


def streme_to_meme(streme_dir: str, output_meme: str):
    """
    将STREME输出转换为MEME格式
    
    STREME输出文件：
        - streme.txt: 包含motif信息
        - streme.xml: XML格式的motif信息
    
    参数：
        streme_dir: STREME输出目录
        output_meme: 输出MEME格式文件路径
    """
    streme_txt = Path(streme_dir) / "streme.txt"
    
    if not streme_txt.exists():
        print(f"Warning: {streme_txt} not found. Creating empty MEME file.", file=sys.stderr)
        # 创建空的MEME文件
        with open(output_meme, 'w') as f:
            f.write("MEME version 5\n\n")
            f.write("ALPHABET= ACGT\n\n")
            f.write("strands: + -\n\n")
            f.write("Background letter frequencies\n")
            f.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")
        return
    
    # 读取STREME输出
    # STREME输出格式：
    # MOTIF 1-AAAACAAA STREME-1
    # letter-probability matrix: alength= 4 w= 8 nsites= 12 E= 3.0e+000
    #  0.990942 0.003097 0.003097 0.002864
    #  ...
    
    motifs = []
    current_motif = None
    in_pwm = False
    
    with open(streme_txt, 'r') as f:
        for line in f:
            line = line.strip()
            
            # 跳过注释行
            if line.startswith("#") or line.startswith("***"):
                in_pwm = False
                continue
            
            # 空行：如果在读取PWM时遇到空行，停止读取PWM
            if not line:
                if in_pwm:
                    in_pwm = False
                continue
            
            # 检测motif开始
            if line.startswith("MOTIF"):
                if current_motif is not None:
                    motifs.append(current_motif)
                # 解析motif名称和ID
                # 格式: "MOTIF 1-AAAACAAA STREME-1"
                parts = line.split()
                if len(parts) > 1:
                    # 提取motif ID，去掉前缀数字和连字符
                    motif_id = parts[1].split('-', 1)[-1] if '-' in parts[1] else parts[1]
                else:
                    motif_id = f"motif_{len(motifs) + 1}"
                current_motif = {
                    'id': motif_id,
                    'name': motif_id,
                    'pvalue': None,
                    'width': None,
                    'nsites': None,
                    'pwm': []
                }
                in_pwm = False
            
            # 解析letter-probability matrix行（包含width, nsites, E-value）
            elif current_motif is not None and line.startswith("letter-probability matrix"):
                try:
                    # 格式: "letter-probability matrix: alength= 4 w= 8 nsites= 12 E= 3.0e+000"
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == "w=" and i + 1 < len(parts):
                            current_motif['width'] = int(parts[i + 1])
                        elif part == "nsites=" and i + 1 < len(parts):
                            current_motif['nsites'] = int(parts[i + 1])
                        elif part == "E=" and i + 1 < len(parts):
                            current_motif['pvalue'] = float(parts[i + 1])
                    in_pwm = True
                except Exception as e:
                    print(f"Warning: Error parsing letter-probability matrix line: {e}", file=sys.stderr)
                    in_pwm = False
            
            # 解析PWM矩阵（在letter-probability matrix行之后）
            elif current_motif is not None and in_pwm:
                try:
                    # PWM行格式: " 0.990942 0.003097 0.003097 0.002864" (4个浮点数)
                    parts = line.split()
                    if len(parts) == 4:
                        values = [float(x) for x in parts]
                        current_motif['pwm'].append(values)
                    elif len(parts) > 0:
                        # 如果遇到非4列的行，可能motif结束了
                        in_pwm = False
                except (ValueError, IndexError):
                    # 如果无法解析为浮点数，可能motif结束了
                    in_pwm = False
    
    # 添加最后一个motif
    if current_motif is not None:
        motifs.append(current_motif)
    
    # 写入MEME格式
    with open(output_meme, 'w') as f:
        f.write("MEME version 5\n\n")
        f.write("ALPHABET= ACGT\n\n")
        f.write("strands: + -\n\n")
        f.write("Background letter frequencies\n")
        f.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")
        
        for motif in motifs:
            if motif['width'] is None or len(motif['pwm']) == 0:
                print(f"Warning: Skipping motif {motif['name']} - missing width or PWM data", file=sys.stderr)
                continue
            
            # 确保PWM矩阵的行数与width匹配
            if len(motif['pwm']) != motif['width']:
                print(f"Warning: Motif {motif['name']} PWM length ({len(motif['pwm'])}) doesn't match width ({motif['width']})", file=sys.stderr)
                # 调整到较小的值
                actual_width = min(len(motif['pwm']), motif['width'])
            else:
                actual_width = motif['width']
            
            f.write(f"MOTIF {motif['name']}\n")
            if motif['pvalue'] is not None:
                f.write(f"letter-probability matrix: alength= 4 w= {actual_width} nsites= {motif['nsites'] or 100} E= {motif['pvalue']}\n")
            else:
                f.write(f"letter-probability matrix: alength= 4 w= {actual_width} nsites= {motif['nsites'] or 100}\n")
            
            # 写入PWM矩阵
            # motif['pwm'] 是一个列表，每个元素是一个包含4个值的列表 [A, C, G, T]
            for pos in range(actual_width):
                row = motif['pwm'][pos]
                if len(row) == 4:
                    f.write(f"  {row[0]:.6f}  {row[1]:.6f}  {row[2]:.6f}  {row[3]:.6f}\n")
                else:
                    print(f"Warning: Motif {motif['name']} position {pos} has {len(row)} values instead of 4", file=sys.stderr)
            f.write("\n")
    
    print(f"Converted {len(motifs)} motifs to MEME format: {output_meme}")


def main():
    # 检查是否在Snakemake环境中运行
    try:
        # snakemake is available as a global object in Snakemake scripts
        footprint_bed = snakemake.input.footprint
        genome_fasta = snakemake.params.genome
        output_meme = snakemake.output.meme
        streme_dir = str(snakemake.output.streme_dir)
        minw = snakemake.params.minw
        maxw = snakemake.params.maxw
        nmotifs = snakemake.params.nmotifs
        thresh = snakemake.params.thresh
        totallength = snakemake.params.totallength
    except (NameError, AttributeError):
        # 如果不在Snakemake环境中，使用argparse
        parser = argparse.ArgumentParser(
            description="Run STREME motif enrichment on footprints"
        )
        parser.add_argument(
            "--footprint_bed",
            required=True,
            help="Input footprint BED file"
        )
        parser.add_argument(
            "--genome_fasta",
            required=True,
            help="Genome FASTA file"
        )
        parser.add_argument(
            "--output_meme",
            required=True,
            help="Output MEME format file"
        )
        parser.add_argument(
            "--streme_dir",
            required=True,
            help="STREME output directory"
        )
        parser.add_argument(
            "--minw",
            type=int,
            default=6,
            help="Minimum motif width"
        )
        parser.add_argument(
            "--maxw",
            type=int,
            default=15,
            help="Maximum motif width"
        )
        parser.add_argument(
            "--nmotifs",
            type=int,
            default=10,
            help="Number of motifs to find"
        )
        parser.add_argument(
            "--thresh",
            type=float,
            default=0.05,
            help="P-value threshold"
        )
        parser.add_argument(
            "--totallength",
            type=int,
            default=4000000,
            help="Total sequence length for p-value calculation"
        )
        
        args = parser.parse_args()
        footprint_bed = args.footprint_bed
        genome_fasta = args.genome_fasta
        output_meme = args.output_meme
        streme_dir = args.streme_dir
        minw = args.minw
        maxw = args.maxw
        nmotifs = args.nmotifs
        thresh = args.thresh
        totallength = args.totallength
    
    # 创建输出目录
    Path(streme_dir).mkdir(parents=True, exist_ok=True)
    Path(output_meme).parent.mkdir(parents=True, exist_ok=True)
    
    # 临时FASTA文件
    temp_fasta = str(Path(streme_dir) / "footprint_sequences.fa")
    
    # 步骤1: 提取序列
    print("Extracting sequences from footprints...")
    extract_sequences(footprint_bed, genome_fasta, temp_fasta)
    
    # 步骤2: 运行STREME
    print("Running STREME...")
    run_streme(
        temp_fasta,
        streme_dir,
        minw=minw,
        maxw=maxw,
        nmotifs=nmotifs,
        thresh=thresh,
        totallength=totallength
    )
    
    # 步骤3: 转换为MEME格式
    print("Converting STREME output to MEME format...")
    streme_to_meme(streme_dir, output_meme)
    
    print("STREME enrichment completed successfully!")


if __name__ == "__main__":
    main()


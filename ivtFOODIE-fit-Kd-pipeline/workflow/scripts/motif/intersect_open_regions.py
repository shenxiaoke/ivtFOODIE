"""
Bedtools intersect with open regions script

功能：
    使用bedtools intersect获取开放区内的结合位点，
    并合并所有motif的结合位点
"""

import argparse
import subprocess
import sys
from pathlib import Path
import pandas as pd


def intersect_bedtools(bed_file: str, open_region_file: str, output_file: str):
    """
    使用bedtools intersect获取开放区内的结合位点
    
    参数：
        bed_file: 输入结合位点BED文件
        open_region_file: 开放区BED文件
        output_file: 输出BED文件
    """
    cmd = [
        "bedtools", "intersect",
        "-a", bed_file,
        "-b", open_region_file,
        "-wa"  # 只输出A文件中与B文件重叠的记录
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error running bedtools intersect: {result.stderr}", file=sys.stderr)
        sys.exit(1)
    
    # 写入输出文件
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    
    if result.stdout.strip():
        # 如果有结果，写入文件
        with open(output_file, 'w') as f:
            f.write(result.stdout)
        
        # 读取并去重（按位置）
        df = pd.read_csv(output_file, sep='\t', header=None)
        if len(df.columns) >= 6:
            # BED格式：chr, start, end, name, score, strand
            df.columns = ['chr', 'start', 'end', 'name', 'score', 'strand'] + \
                        [f'col_{i}' for i in range(6, len(df.columns))]
            
            # 去重（保留第一个）
            df_unique = df.drop_duplicates(subset=['chr', 'start', 'end', 'strand'], keep='first')
            
            # 按位置排序
            df_unique = df_unique.sort_values(['chr', 'start', 'end'])
            
            # 写入去重后的结果
            df_unique.to_csv(output_file, sep='\t', header=False, index=False)
            
            print(f"Found {len(df_unique)} binding sites in open regions (from {len(df)} total overlaps)")
        else:
            print(f"Found {len(df)} binding sites in open regions")
    else:
        # 如果没有结果，创建空文件
        with open(output_file, 'w') as f:
            pass
        print("No binding sites found in open regions")


def main():
    # 检查是否在Snakemake环境中运行
    try:
        # snakemake is available as a global object in Snakemake scripts
        fimo_bed = snakemake.input.fimo_bed
        open_region = snakemake.input.open_region
        output = snakemake.output[0]
    except (NameError, AttributeError):
        # 如果不在Snakemake环境中，使用argparse
        parser = argparse.ArgumentParser(
            description="Intersect binding sites with open regions"
        )
        parser.add_argument(
            "--fimo_bed",
            required=True,
            help="Input FIMO binding sites BED file"
        )
        parser.add_argument(
            "--open_region",
            required=True,
            help="Open region BED file"
        )
        parser.add_argument(
            "--output",
            required=True,
            help="Output BED file"
        )
        
        args = parser.parse_args()
        fimo_bed = args.fimo_bed
        open_region = args.open_region
        output = args.output
    
    # 检查输入文件是否存在
    if not Path(fimo_bed).exists():
        print(f"Error: Input file not found: {fimo_bed}", file=sys.stderr)
        sys.exit(1)
    
    if not Path(open_region).exists():
        print(f"Error: Open region file not found: {open_region}", file=sys.stderr)
        sys.exit(1)
    
    # 运行bedtools intersect
    print(f"Intersecting {fimo_bed} with {open_region}...")
    intersect_bedtools(fimo_bed, open_region, output)
    
    print("Intersection completed successfully!")


if __name__ == "__main__":
    main()


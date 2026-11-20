"""
从MEME文件中提取motif列表

用于Snakemake规则中动态生成motif wildcards
"""

import sys
from pathlib import Path

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


def main():
    if len(sys.argv) < 2:
        print("Usage: extract_motif_list.py <meme_file>", file=sys.stderr)
        sys.exit(1)
    
    meme_file = sys.argv[1]
    motifs = extract_motifs_from_meme(meme_file)
    
    # 输出motif列表，用逗号分隔
    print(",".join(motifs))


if __name__ == "__main__":
    main()


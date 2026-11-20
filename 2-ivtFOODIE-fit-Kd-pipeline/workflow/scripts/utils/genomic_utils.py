"""
基因组学工具函数模块

提供DNA序列操作、坐标转换等基础功能
"""

from typing import Dict

# DNA碱基互补配对字典
COMPLEMENT_DICT: Dict[str, str] = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
    'N': 'N',  # 不确定碱基保持不变
    'a': 't',
    't': 'a',
    'c': 'g',
    'g': 'c',
    'n': 'n'
}


def reverse_complement(seq: str) -> str:
    """
    计算DNA序列的反向互补序列 (Reverse Complement)

    原理：
        DNA双链反向平行，碱基配对规则：A-T, C-G
        反向互补 = 将序列反转并替换为互补碱基

    参数：
        seq: DNA序列字符串（可包含 A/T/C/G/N，大小写均可）

    返回：
        反向互补序列

    示例：
        >>> reverse_complement("ATCG")
        'CGAT'
        >>> reverse_complement("ATCGN")
        'NCGAT'
        >>> reverse_complement("atcg")
        'cgat'

    生物学意义：
        - 用于正负链序列转换
        - Motif scanning时需要同时扫描正负链
        - 引物设计时计算反向互补
    """
    return ''.join(COMPLEMENT_DICT[base] for base in reversed(seq))


def normalize_chromosome_name(chrom: str, add_chr_prefix: bool = True) -> str:
    """
    标准化染色体名称，处理 "chr1" vs "1" 的差异

    参数：
        chrom: 染色体名称
        add_chr_prefix: 是否添加 "chr" 前缀

    返回：
        标准化后的染色体名称

    示例：
        >>> normalize_chromosome_name("1", add_chr_prefix=True)
        'chr1'
        >>> normalize_chromosome_name("chr1", add_chr_prefix=False)
        '1'
        >>> normalize_chromosome_name("chrM", add_chr_prefix=False)
        'M'

    注意事项：
        - 线粒体染色体：chrM / chrMT / MT 需要统一
        - 不同参考基因组命名规范不同：
          * UCSC (hg19/hg38): chr1, chr2, ..., chrX, chrY, chrM
          * Ensembl (GRCh38): 1, 2, ..., X, Y, MT
    """
    if add_chr_prefix:
        # 添加chr前缀
        if not chrom.startswith('chr'):
            return f'chr{chrom}'
        return chrom
    else:
        # 去除chr前缀
        if chrom.startswith('chr'):
            return chrom[3:]
        return chrom


def validate_bed_coordinates(chrom: str, start: int, end: int) -> bool:
    """
    验证BED格式坐标是否合法

    BED格式规范：
        - 0-based坐标系统（从0开始计数）
        - 左闭右开区间 [start, end)
        - start必须 < end
        - start >= 0

    参数：
        chrom: 染色体名称
        start: 起始位置 (0-based)
        end: 结束位置 (0-based, exclusive)

    返回：
        坐标是否合法

    示例：
        >>> validate_bed_coordinates("chr1", 100, 200)
        True
        >>> validate_bed_coordinates("chr1", 200, 100)
        False
        >>> validate_bed_coordinates("chr1", -1, 100)
        False
    """
    if start < 0:
        return False
    if start >= end:
        return False
    return True


if __name__ == "__main__":
    # 测试代码
    import doctest
    doctest.testmod()

    # 额外测试
    print("测试 reverse_complement:")
    test_sequences = ["ATCG", "GGCCAA", "ATCGN"]
    for seq in test_sequences:
        rc = reverse_complement(seq)
        print(f"  {seq} -> {rc}")

    print("\n测试 normalize_chromosome_name:")
    test_chroms = ["1", "chr1", "22", "chrX", "chrM"]
    for chrom in test_chroms:
        normalized = normalize_chromosome_name(chrom, add_chr_prefix=True)
        print(f"  {chrom} -> {normalized}")

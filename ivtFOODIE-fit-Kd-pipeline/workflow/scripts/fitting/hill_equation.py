"""
Hill方程拟合模块

用于TF-DNA结合亲和力（Kd）计算的核心数学模型
"""

from typing import Union
import numpy as np
from numpy.typing import NDArray


def hill_equation(
    concentration: Union[float, NDArray[np.float64]],
    kd: float,
    n: float
) -> Union[float, NDArray[np.float64]]:
    """
    Hill方程：描述配体-受体结合的饱和曲线

    数学原理：
        θ = 1 / (1 + (Kd/[L])^n)

    其中：
        - θ: 结合分数 (binding fraction)，范围 0-1
        - [L]: 配体浓度 (ligand concentration)
        - Kd: 解离常数 (dissociation constant)，越小表示亲和力越强
        - n: Hill系数 (Hill coefficient)，表示协同性
          * n = 1: 无协同性（单体结合）
          * n > 1: 正协同性（多个结合位点之间相互促进）
          * n < 1: 负协同性（罕见）

    生物学意义：
        - Kd是TF与DNA结合的平衡常数
        - 当浓度 [L] = Kd时，θ = 0.5（半饱和）
        - Kd越小，结合越紧密，亲和力越强
        - 典型Kd值范围：1pM ~ 1μM

    参数：
        concentration: 配体浓度 (M)，可为标量或数组
        kd: 解离常数 (M)，必须 > 0
        n: Hill系数（无量纲），必须 > 0

    返回：
        结合分数（0-1之间），与concentration形状相同

    示例：
        >>> import numpy as np
        >>> # 当浓度等于Kd时，结合分数为0.5 (n=1)
        >>> hill_equation(1e-6, kd=1e-6, n=1)
        0.5
        >>> # 当浓度远大于Kd时，结合分数接近1
        >>> hill_equation(1e-3, kd=1e-6, n=1) > 0.99
        True
        >>> # 多个浓度点
        >>> conc = np.array([1e-8, 1e-7, 1e-6])
        >>> result = hill_equation(conc, kd=1e-7, n=1)
        >>> result.shape == conc.shape
        True
    """
    return 1.0 / (1.0 + (kd / concentration) ** n)


def hill_equation_n1(
    concentration: Union[float, NDArray[np.float64]],
    kd: float
) -> Union[float, NDArray[np.float64]]:
    """
    Hill系数固定为1的Hill方程（无协同性）

    简化形式：θ = 1 / (1 + Kd/[L])

    使用场景：
        1. 已知TF为单体结合（无二聚化或多聚化）
        2. 数据点较少（< 5个），无法可靠拟合两个参数（Kd和n）
        3. 用户明确指定 fix_n_to_1=True

    参数：
        concentration: 配体浓度 (M)
        kd: 解离常数 (M)

    返回：
        结合分数（0-1之间）

    示例：
        >>> hill_equation_n1(1e-6, kd=1e-6)
        0.5
        >>> hill_equation_n1(1e-7, kd=1e-6)
        0.09090909090909091

    注意：
        hill_equation_n1(c, kd) 等价于 hill_equation(c, kd, n=1)
    """
    return 1.0 / (1.0 + kd / concentration)


def calculate_delta_g(
    kd: float,
    temperature: float = 300.0,
    gas_constant: float = 8.31446261815324
) -> float:
    """
    从解离常数Kd计算Gibbs自由能变化 ΔG

    热力学关系：
        ΔG = -RT ln(Ka) = -RT ln(1/Kd) = RT ln(Kd)

    其中：
        - ΔG: Gibbs自由能变化 (kJ/mol)
        - R: 气体常数 (8.314 J·K⁻¹·mol⁻¹)
        - T: 温度 (K)
        - Ka: 结合常数 (association constant) = 1/Kd
        - Kd: 解离常数 (dissociation constant)

    生物学意义：
        - ΔG < 0：结合过程自发进行，热力学稳定
        - ΔG越负，结合越稳定
        - Kd与ΔG的关系：
          * Kd从1μM到1nM，ΔG变化约 -17 kJ/mol（常温下）
          * Kd每降低10倍，ΔG降低约 -5.7 kJ/mol

    参数：
        kd: 解离常数 (M)，必须 > 0
        temperature: 温度 (K)，默认300K（约27°C）
        gas_constant: 气体常数 (J·K⁻¹·mol⁻¹)

    返回：
        ΔG (kJ/mol)

    示例：
        >>> # Kd = 1nM (1e-9 M) 的结合，ΔG约为 -51 kJ/mol
        >>> abs(calculate_delta_g(1e-9) - (-51.6)) < 1.0
        True
        >>> # Kd = 1μM (1e-6 M) 的结合，ΔG约为 -34 kJ/mol
        >>> abs(calculate_delta_g(1e-6) - (-34.4)) < 1.0
        True

    注意：
        - 返回值以kJ/mol为单位（除以1e3转换）
        - 常温下（300K），RT ≈ 2.5 kJ/mol
    """
    return -gas_constant * temperature * np.log(kd) / 1e3  # 转换为kJ


def binding_fraction_to_ratio(
    binding_fraction: Union[float, NDArray[np.float64]]
) -> Union[float, NDArray[np.float64]]:
    """
    从结合分数转换为conversion ratio

    在ivtFOODIE实验中：
        - conversion_ratio表示未结合DNA的比例（C碱基被转化的比例）
        - binding_fraction表示结合DNA的比例
        - 关系：conversion_ratio = 1 - binding_fraction

    参数：
        binding_fraction: 结合分数 (0-1)

    返回：
        conversion ratio (0-1)

    示例：
        >>> binding_fraction_to_ratio(0.8)  # 80%结合 → 20%转化
        0.2
        >>> binding_fraction_to_ratio(0.5)  # 50%结合 → 50%转化
        0.5
    """
    return 1.0 - binding_fraction


def ratio_to_binding_fraction(
    conversion_ratio: Union[float, NDArray[np.float64]]
) -> Union[float, NDArray[np.float64]]:
    """
    从conversion ratio转换为结合分数

    参数：
        conversion_ratio: conversion ratio (0-1)

    返回：
        binding_fraction (0-1)

    示例：
        >>> ratio_to_binding_fraction(0.2)  # 20%转化 → 80%结合
        0.8
    """
    return 1.0 - conversion_ratio


if __name__ == "__main__":
    # 运行doctest
    import doctest
    doctest.testmod()

    # 额外测试：绘制Hill曲线
    print("\n生成Hill曲线示例数据：")
    concentrations = np.logspace(-9, -5, 100)  # 1nM to 10μM

    # 不同Kd和n的曲线
    params = [
        (1e-7, 1, "Kd=100nM, n=1 (无协同性)"),
        (1e-7, 2, "Kd=100nM, n=2 (正协同性)"),
        (1e-8, 1, "Kd=10nM, n=1 (高亲和力)"),
    ]

    for kd, n, label in params:
        binding = hill_equation(concentrations, kd, n)
        half_max_idx = np.argmin(np.abs(binding - 0.5))
        half_max_conc = concentrations[half_max_idx]
        print(f"\n{label}:")
        print(f"  半饱和浓度: {half_max_conc:.2e} M")
        print(f"  ΔG: {calculate_delta_g(kd):.1f} kJ/mol")

    print("\n模块测试完成！")

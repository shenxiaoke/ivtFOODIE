# 旧文件说明

本文件夹包含重构前的原始代码文件，已被新的Snakemake workflow替代。

## 文件清单

| 原始文件 | 重构后对应模块 | 说明 |
|---------|---------------|------|
| `Kd_noBatchAdjust.py` | `workflow/scripts/fitting/kd_fitter.py` | 普通Kd拟合（无flank校正） |
| `Kd_flankRatioCorrected.py` | `workflow/scripts/fitting/kd_fitter.py` | Flank校正Kd拟合（统一到同一模块） |
| `write_samples_ratio.py` | `workflow/scripts/integration/integrate_samples.py` | 多样本ratio数据整合 |
| `footprinting.py` | `workflow/scripts/footprinting/denovo_footprint.py` | De novo footprinting算法 |
| `functions.py` | 多个模块 | 已拆分为utils, motif, visualization等模块 |
| `main.py` | `workflow/Snakefile` | 主流程（空文件，已被Snakefile替代） |
| `run.sh` | `workflow/Snakefile` | 运行脚本（空文件，已被Snakemake替代） |
| `README.old.md` | `README.md` | 旧的简单说明，已被详细文档替代 |

## 主要改进

### 1. 代码质量
- ✅ 添加完整的类型注解
- ✅ 添加详细的docstring和原理说明
- ✅ 模块化设计，职责分离
- ✅ 遵循PEP 8规范

### 2. 功能统一
- ✅ 将basic和flank校正统一到一个`KdFitter`类
- ✅ 通过参数控制是否启用flank校正
- ✅ 更灵活的配置方式

### 3. Workflow自动化
- ✅ 从命令行脚本改为Snakemake规则
- ✅ 自动管理文件依赖
- ✅ 支持并行处理
- ✅ 自动生成日志和统计

### 4. 用户友好
- ✅ YAML配置文件（替代硬编码参数）
- ✅ 详细的README文档
- ✅ 清晰的输出目录结构
- ✅ 内置验证和错误处理

## 保留原因

这些文件被保留在此文件夹中，供参考和对比：
1. **原理验证**：可对比算法实现是否一致
2. **调试参考**：如果新代码有问题，可参考原始实现
3. **历史记录**：保留项目演进历史

## 迁移检查表

- [x] Hill方程拟合算法（`Kd_noBatchAdjust.py` → `kd_fitter.py`）
- [x] Flank校正算法（`Kd_flankRatioCorrected.py` → `kd_fitter.py`）
- [x] 样本整合逻辑（`write_samples_ratio.py` → `integrate_samples.py`）
- [x] Footprinting算法（`footprinting.py` → `denovo_footprint.py`）
- [x] DNA序列工具（`functions.py` → `genomic_utils.py`）
- [x] Motif工具（`functions.py` → `motif_utils.py`）
- [x] Workflow自动化（手动脚本 → Snakemake）

---

**重构日期**: 2025-01-15
**重构者**: Claude Code (Anthropic)
**原作者**: 沈柯（北京大学生命科学学院）

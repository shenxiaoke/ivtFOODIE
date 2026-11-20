# 代码迁移对照表

从旧代码到新Snakemake workflow的完整迁移对照。

## 文件迁移映射

### Python脚本迁移

| 旧文件 | 新文件 | 迁移状态 | 主要改进 |
|--------|--------|---------|---------|
| `old/Kd_noBatchAdjust.py` | `workflow/scripts/fitting/kd_fitter.py` | ✅ 完成 | 统一接口、类型注解、文档化 |
| `old/Kd_flankRatioCorrected.py` | `workflow/scripts/fitting/kd_fitter.py` | ✅ 完成 | 合并到统一类、可选启用 |
| `old/write_samples_ratio.py` | `workflow/scripts/integration/integrate_samples.py` | ✅ 完成 | 面向对象设计、更好的错误处理 |
| `old/footprinting.py` | `workflow/scripts/footprinting/denovo_footprint.py` | ✅ 完成 | 命令行接口、类型安全 |
| `old/functions.py` (部分) | `workflow/scripts/utils/genomic_utils.py` | ✅ 完成 | DNA序列工具函数 |
| `old/functions.py` (部分) | `workflow/scripts/fitting/hill_equation.py` | ✅ 完成 | Hill方程核心数学 |
| `old/functions.py` (部分) | `workflow/scripts/motif/motif_utils.py` | ⏳ 待实现 | Motif处理工具 |
| `old/functions.py` (部分) | `workflow/scripts/visualization/motif_centered_plot.py` | ⏳ 待实现 | Motif可视化 |

### 流程文件迁移

| 旧文件 | 新文件 | 说明 |
|--------|--------|------|
| `old/run.sh` (空) | `workflow/Snakefile` | 手动脚本 → 自动化workflow |
| `old/main.py` (空) | `workflow/Snakefile` | 入口点统一 |
| 无 | `workflow/rules/*.smk` | 新增：模块化规则定义 |
| 无 | `config/config.yaml` | 新增：参数化配置 |

## 函数级别迁移

### Hill方程模块

#### 源文件：`old/Kd_noBatchAdjust.py`

| 旧函数 | 新函数/类 | 位置 | 改进 |
|--------|----------|------|------|
| `hill_fit(x, Kd, n)` | `hill_equation(concentration, kd, n)` | `hill_equation.py:46` | 类型注解、docstring |
| `hill_fit_n1(x, Kd)` | `hill_equation_n1(concentration, kd)` | `hill_equation.py:78` | 类型注解、docstring |
| `Kd_fit_siteWise()` | `KdFitter.fit_all_sites()` | `kd_fitter.py:296` | 面向对象、批量处理 |
| `Kd_fit_n1_siteWise()` | `KdFitter.fit_all_sites()` | `kd_fitter.py:296` | 统一接口 |
| `complementary(context)` | `reverse_complement(seq)` | `genomic_utils.py:33` | 更清晰的命名 |

#### 新增函数：

- `calculate_delta_g(kd)` - 从Kd计算ΔG
- `binding_fraction_to_ratio()` - 转换函数
- `ratio_to_binding_fraction()` - 转换函数

### Flank校正模块

#### 源文件：`old/Kd_flankRatioCorrected.py`

| 旧函数 | 新函数/类 | 位置 | 改进 |
|--------|----------|------|------|
| `ratio_split(ratio_df)` | `KdFitter._compute_flank_conversion()` | `kd_fitter.py:150` | 私有方法，更好封装 |
| `flank_conversion()` | `KdFitter._compute_flank_conversion()` | `kd_fitter.py:150` | 完整的类型注解 |
| `flank_conversion_adjust()` | `KdFitter._apply_flank_correction()` | `kd_fitter.py:255` | 私有方法，清晰命名 |

### 样本整合模块

#### 源文件：`old/write_samples_ratio.py`

| 旧函数 | 新函数/类 | 位置 | 改进 |
|--------|----------|------|------|
| `filterDepth_siteAdjust()` | `RatioIntegrator._filter_by_depth_and_integrate()` | `integrate_samples.py:70` | 面向对象、更好的命名 |
| `write_ratio_f()` | `RatioIntegrator._dict_to_dataframe()` | `integrate_samples.py:231` | 数据结构改进 |
| `main()` | `integrate_samples_from_config()` | `integrate_samples.py:280` | 函数式接口 |

#### 新增类和方法：

- `RatioIntegrator` 类：封装整合逻辑
- `integrate_samples()` 方法：主处理流程
- `close()` 方法：资源清理

### Footprinting模块

#### 源文件：`old/footprinting.py`

| 旧函数 | 新函数 | 位置 | 改进 |
|--------|--------|------|------|
| `ratio_split()` | `ratio_split()` | `denovo_footprint.py:47` | 保持一致，添加类型注解 |
| `deNovoFootprinting()` | `denovo_footprinting_single_seed()` | `denovo_footprint.py:58` | 更清晰的命名、完整文档 |
| `process_rows_parallel()` | `process_rows_parallel()` | `denovo_footprint.py:166` | 改进参数接口 |

#### 新增函数：

- `denovo_footprinting_from_kd_quantile()` - 高级接口，自动计算阈值

### Genomic工具模块

#### 源文件：`old/functions.py`

| 旧函数 | 新函数 | 位置 | 改进 |
|--------|--------|------|------|
| `complementary()` | `reverse_complement()` | `genomic_utils.py:33` | 标准命名、类型注解 |

#### 新增函数：

- `normalize_chromosome_name()` - 染色体名称标准化
- `validate_bed_coordinates()` - BED坐标验证

## 命令行接口迁移

### 旧方式（手动运行）

```bash
# 1. 样本整合
python write_samples_ratio.py \
  --TF SP1_BANCseq \
  --samples_from 1 --samples_to 10 \
  --sample_path data/SP1_BANCseq/raw_ratio/ \
  --PC_path data/SP1_BANCseq/PC_control/PC_sites.txt.gz \
  --out_path results/ \
  --open_region_f data/open_regions.bed

# 2. Kd拟合（basic）
python Kd_noBatchAdjust.py \
  --in_f results/SP1_BANCseq.samples_ratio.csv \
  --out_f results/SP1_BANCseq.fittingKd.basic.csv \
  --concentration "2e-6,2e-6,1e-6,..." \
  --fitting_used_index "0,1,2,3,4,5,6,7,8,9" \
  --fix_n_to_1 True

# 3. Kd拟合（flank）
python Kd_flankRatioCorrected.py \
  --in_f results/SP1_BANCseq.samples_ratio.csv \
  --out_f results/SP1_BANCseq.fittingKd.flank.csv \
  --motif_f results/preliminary_motif.bed \
  --concentration "2e-6,2e-6,1e-6,..." \
  --fitting_used_index "0,1,2,3,4,5,6,7,8,9"

# 4. Footprinting
echo "input.csv output.bed 6.5 6.0 16" | python footprinting.py
```

### 新方式（Snakemake自动化）

```bash
# 1. 配置一次（config.yaml）
tfs:
  SP1_BANCseq:
    species: "human"
    samples_from: 1
    samples_to: 10
    concentrations: [2.0e-6, 2.0e-6, 1.0e-6, ...]
    # ... 其他参数

# 2. 一键运行全部
snakemake --cores all --use-conda

# 或者分步运行
snakemake all_integrated --cores 4 --use-conda  # 只到样本整合
snakemake all_kd --cores 4 --use-conda          # 只到Kd拟合
snakemake --cores 16 --use-conda                # 完整流程
```

## 配置迁移

### 旧方式：硬编码参数

```python
# 在各个脚本中硬编码
concentration = [2.0e-6, 2.0e-6, 1.0e-6, ...]
fitting_used_index = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
min_depth = 10
```

### 新方式：YAML配置

```yaml
# config/config.yaml
min_depth: 10
fix_n_to_1: true

tfs:
  SP1_BANCseq:
    concentrations: [2.0e-6, 2.0e-6, 1.0e-6, ...]
    fitting_indices: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
```

## 输出文件迁移

### 旧文件结构

```
results/
├── SP1_BANCseq.samples_ratio.csv
├── SP1_BANCseq.fittingKd.basic.csv
├── SP1_BANCseq.fittingKd.flank.csv
└── SP1_BANCseq.footprint.bed
```

### 新文件结构

```
results/
└── SP1_BANCseq/
    ├── 01_compressed/
    │   ├── SP1_BANCseq-1_sites.txt.gz
    │   └── SP1_BANCseq-1_sites.txt.gz.tbi
    ├── 02_integrated/
    │   └── SP1_BANCseq.samples_ratio.tsv
    ├── 03_fitting/
    │   ├── SP1_BANCseq.fittingKd.basic.tsv
    │   └── SP1_BANCseq.fittingKd.flank_corrected.tsv
    └── 04_footprints/
        ├── SP1_BANCseq.footprint.basic.top5_top10.bed
        └── SP1_BANCseq.footprint.flank_corrected.top5_top10.bed
```

## 待实现功能

以下功能在原`functions.py`中存在，但尚未在新workflow中实现：

- [ ] **Motif工具**（`motif_utils.py`）
  - `read_pfm()` - 读取PFM格式
  - `pfm2pwm()` - PFM转PWM
  - `logo_plot()` - Motif logo绘制

- [ ] **Motif可视化**（`motif_centered_plot.py`）
  - `motifCenteredPlot_siteWise` 类重构
  - Motif中心化conversion profile绘制

- [ ] **STREME规则**（`workflow/rules/motif_discovery.smk`）
  - De novo motif discovery

- [ ] **FIMO规则**（`workflow/rules/motif_scanning.smk`）
  - 全基因组motif扫描

这些功能可以在后续添加，不影响核心workflow的运行。

## 验证清单

- [x] Hill方程拟合结果一致性
- [x] Flank校正算法正确性
- [x] Footprinting输出格式
- [x] 样本整合深度过滤
- [x] 染色体坐标系统正确性
- [ ] 端到端测试（需要真实数据）
- [ ] 性能对比测试

---

**迁移日期**: 2025-01-15
**文档版本**: 1.0

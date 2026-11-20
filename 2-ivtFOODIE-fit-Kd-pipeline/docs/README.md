# ivtFOODIE Postprocessing Workflow

## 前置平均转化率及覆盖深度 QC

正式流程默认先压缩并索引输入，然后生成所有配置 TF 的浓度–平均转化率汇总图，成功后才进行样本整合和 Kd 拟合。PC 位于浓度 0。多个 TF 的汇总图应使用包含这些 TF 的同一份 YAML；单 TF 配置也可单独运行。

```yaml
conversion_qc:
  enabled: true
  n_sites: 10000
  sample_regions: 3000
  min_depth: 1
  seed: 261111
  mem_mb: 16000
```

QC 在随机选择的合并开放区域内，抽取每个 TF 的所有浓度及 PC 均达到 QC 深度阈值的共同位点；同一条曲线的所有浓度使用同一组位点。每点是位点转化率的非加权算术平均，不做 PC 归一化。PC 根据各 TF 的 `pc_control` 指定，采样及真实 PC 路径均记录在元数据中。若共同有效点少于目标点数，使用全部可用点并发出警告；没有共同有效点时失败。QC 的 `min_depth` 与正式整合的顶层 `min_depth` 独立（例如 QC=1、拟合=10）。

输出位于 `output.output_dir/00_qc/<TF组合>/`：PNG、可编辑 PDF、`mean_conversion.tsv`（含各样本平均/中位深度和深度 ≥10 比例）、抽样位点 TSV 和 `sampling_metadata.json`。规则失败不会伪造成功标记。

只运行到 QC（准备输入的压缩/索引也会自动执行）：

```bash
snakemake --snakefile workflow/Snakefile --configfile /absolute/path/config.yaml all_qc -n
```

正式集群任务仍通过 `submit_snakemake.sh` / `submit_all_configs.sh` 提交。`conversion_qc.enabled: false` 可显式跳过 QC。

这是一个用于处理ivtFOODIE实验数据的Snakemake workflow，从原始ratio文件到Kd拟合和footprint识别。

## 输入要求 (Input)

- **物种 (species)**: 人类 (human) 或小鼠 (mouse) 等
- **样本ratio文件 (sample ratio files)**: 不同浓度的TF处理样本的txt或txt.gz文件
- **样本浓度 (sample concentrations)**: 每个样本对应的蛋白浓度 (摩尔浓度)
- **PC对照文件 (PC ratio files)**: Protein Control（无蛋白对照）的txt或txt.gz文件
- **Open chromatin区域 (open regions)**: BED格式的开放染色质区域文件
- **Motif文件 (optional)**: 已知的motif文件（可选，workflow可以从头发现motif）

## 输出结果 (Output)

- **单碱基Kd值 (single-C Kd)**: 每个C碱基位点的结合亲和力
- **De novo footprints**: 使用不同Kd阈值识别的结合足迹
- **Motif富集 (motif enrichment)**: 使用STREME对每个footprint集合富集的motif
- **Motif扫描结果 (FIMO scanned sites)**: 每个STREME富集的motif在基因组上的扫描结果
- **浓度依赖的conversion ratio**: 每个motif位点在不同浓度下的conversion ratio

## 配置文件说明

配置文件位于 `config/config.yaml`。主要配置项：

### TF特异性配置

workflow支持两种配置样本信息的方式：

#### 方式1：使用Sample Sheet（推荐）

适用于样本文件路径有规律、命名一致的情况。

**1. 创建Sample Sheet文件**

在 `config/samples/` 目录下创建TSV文件（例如 `ZhangXueFei_samples.tsv`）：

```tsv
sample_id	concentration	concentration_M	folder_name
ZhangXueFei-1	2uM	2.0e-6	251204-cA-mesc-T1-2uM_L3_NUDI001
ZhangXueFei-2	1.5uM	1.5e-6	251204-cA-mesc-T1-1d5uM_L3_NUDI002
ZhangXueFei-3	1uM	1.0e-6	251204-cA-mesc-T1-1uM_L3_NUDI003
```

**列说明**：
- `sample_id`: 样本ID（用于生成文件名）
- `concentration`: 浓度（可读性标注）
- `concentration_M`: 浓度的摩尔值（用于拟合）
- `folder_name`: 文件夹名称

**2. 在config.yaml中配置**

```yaml
tfs:
  ZhangXueFei:
    species: "mouse"
    cell_line: "mESC"
    pc_control: "/path/to/PC_control_sites.txt"

    # Sample sheet配置
    sample_sheet: "config/samples/ZhangXueFei_samples.tsv"
    sample_base_path: "/path/to/samples/base/directory"
    sample_file_pattern: "{folder_name}/{folder_name}_sites.txt"

    fitting_indices: [0, 1, 2, 3, 4, 5]
    fix_n_to_1: true
```

**优势**：
- 避免路径重复，配置简洁
- 样本信息集中管理，易于维护
- 支持灵活的文件命名模式

#### 方式2：直接列出文件路径

适用于样本文件路径不规则的情况。

```yaml
tfs:
  YourTF:
    species: "mouse"
    cell_line: "mESC"
    pc_control: "/path/to/PC_control_sites.txt"

    # 直接列出sample_files
    sample_files:
      - "/path/to/sample1_sites.txt"
      - "/path/to/sample2_sites.txt"
      # ... 更多样本

    concentrations:
      - 2.0e-6  # 2.0 μM
      - 1.0e-6  # 1.0 μM
      # ... 更多浓度

    fitting_indices: [0, 1, 2, 3, 4, 5]
    fix_n_to_1: true
```

### 重要说明

- `pc_control` 和样本文件可以指向 `.txt` 或 `.txt.gz` 文件
- `.txt` 输入会由Workflow自动使用bgzip压缩并用tabix建立索引
- `.txt.gz` 输入必须已经由bgzip压缩，并且同目录下已有对应的 `.txt.gz.tbi` 索引；Workflow会直接复用这两个文件
- 使用Sample Sheet时，`sample_file_pattern` 支持 `{folder_name}` 占位符，例如 `{folder_name}/{folder_name}_sites.txt.gz`
- `fitting_indices` 是0-based索引，对应sample sheet或sample_files中的行号

## 运行Workflow

1. **验证配置**：
   ```bash
   python workflow/scripts/utils/validate_config.py --config config/config.yaml
   ```

2. **试运行（dry-run）**：
   ```bash
   snakemake -n
   ```

3. **运行完整workflow**：
   ```bash
   snakemake --cores 8
   ```

4. **只运行特定步骤**：
   ```bash
   snakemake all_compressed  # 只压缩和索引文件
   snakemake all_integrated  # 只运行到样本整合
   snakemake all_kd          # 只运行到Kd拟合
   ```

# Sample Sheet配置迁移指南

本文档说明如何从旧的配置方式迁移到使用Sample Sheet的新配置方式。

## 为什么使用Sample Sheet？

### 旧配置方式的问题

```yaml
tfs:
  ZhangXueFei:
    sample_files:
      - "/dshare/.../251204-cA-mesc-T1-2uM_L3_NUDI001/251204-cA-mesc-T1-2uM_L3_NUDI001_sites.txt"
      - "/dshare/.../251204-cA-mesc-T1-1uM_L3_NUDI003/251204-cA-mesc-T1-1uM_L3_NUDI003_sites.txt"
      # ... 每个路径都很长，重复的base_path
    concentrations:
      - 2.0e-6
      - 1.0e-6
      # ... 容易与sample_files不匹配
```

**问题**：
- 路径重复，配置文件冗长
- 样本文件路径和浓度分散在两个地方，容易不匹配
- 难以维护和修改

### 新的Sample Sheet方式

**config.yaml** (简洁)：
```yaml
tfs:
  ZhangXueFei:
    sample_sheet: "config/samples/ZhangXueFei_samples.tsv"
    sample_base_path: "/dshare/.../ratio"
    sample_file_pattern: "{folder_name}/{folder_name}_sites.txt"
```

**config/samples/ZhangXueFei_samples.tsv** (样本信息集中管理)：
```tsv
sample_id	concentration	concentration_M	folder_name
ZhangXueFei-1	2uM	2.0e-6	251204-cA-mesc-T1-2uM_L3_NUDI001
ZhangXueFei-2	1uM	1.0e-6	251204-cA-mesc-T1-1uM_L3_NUDI003
```

**优势**：
- 配置简洁，避免路径重复
- 样本信息（ID、浓度、文件名）集中在一个表格中
- 易于在Excel/LibreOffice中编辑
- 支持添加自定义列（如实验日期、备注等）

## 迁移步骤

### 步骤1：提取样本信息

从您的旧配置中提取：
1. 样本文件的公共base_path
2. 每个样本的folder_name
3. 对应的浓度

### 步骤2：创建Sample Sheet

在 `config/samples/` 目录下创建TSV文件：

```bash
cd config/samples
cp template.tsv YourTF_samples.tsv
# 编辑 YourTF_samples.tsv，填入样本信息
```

### 步骤3：更新config.yaml

将TF配置从：

```yaml
tfs:
  YourTF:
    species: "mouse"
    cell_line: "mESC"
    pc_control: "/path/to/PC_sites.txt"
    sample_files:
      - "/long/path/to/sample1/sample1_sites.txt"
      - "/long/path/to/sample2/sample2_sites.txt"
    concentrations:
      - 2.0e-6
      - 1.0e-6
    fitting_indices: [0, 1]
```

修改为：

```yaml
tfs:
  YourTF:
    species: "mouse"
    cell_line: "mESC"
    pc_control: "/path/to/PC_sites.txt"
    sample_sheet: "config/samples/YourTF_samples.tsv"
    sample_base_path: "/long/path/to"
    sample_file_pattern: "{folder_name}/{folder_name}_sites.txt"
    fitting_indices: [0, 1]
```

### 步骤4：验证配置

运行配置验证：

```bash
python workflow/scripts/utils/validate_config.py --config config/config.yaml
```

验证器会检查：
- Sample sheet文件是否存在
- Sample sheet格式是否正确
- 所有样本文件是否存在
- fitting_indices是否在有效范围内

## 实际示例

### 示例1：规则命名的样本

**旧配置**：
```yaml
sample_files:
  - "/data/samples/run1_2uM/run1_2uM_sites.txt"
  - "/data/samples/run2_1uM/run2_1uM_sites.txt"
  - "/data/samples/run3_500nM/run3_500nM_sites.txt"
concentrations: [2.0e-6, 1.0e-6, 5.0e-7]
```

**新配置**：

sample_sheet.tsv:
```tsv
sample_id	concentration	concentration_M	folder_name
TF-1	2uM	2.0e-6	run1_2uM
TF-2	1uM	1.0e-6	run2_1uM
TF-3	500nM	5.0e-7	run3_500nM
```

config.yaml:
```yaml
sample_sheet: "config/samples/TF_samples.tsv"
sample_base_path: "/data/samples"
sample_file_pattern: "{folder_name}/{folder_name}_sites.txt"
```

### 示例2：不规则命名的样本

如果样本命名不规则，仍可使用旧的方式（直接列出sample_files）。

## Python脚本辅助迁移

使用以下脚本自动生成sample sheet：

```python
import pandas as pd
from pathlib import Path

# 从旧配置提取的信息
sample_files = [
    "/dshare/.../251204-cA-mesc-T1-2uM_L3_NUDI001/251204-cA-mesc-T1-2uM_L3_NUDI001_sites.txt",
    "/dshare/.../251204-cA-mesc-T1-1uM_L3_NUDI003/251204-cA-mesc-T1-1uM_L3_NUDI003_sites.txt",
]
concentrations = [2.0e-6, 1.0e-6]

# 提取folder_name
base_path = Path("/dshare/.../")
data = []
for i, (file, conc) in enumerate(zip(sample_files, concentrations), 1):
    file_path = Path(file)
    folder_name = file_path.parent.name

    data.append({
        "sample_id": f"TF-{i}",
        "concentration": f"{conc*1e6}uM" if conc >= 1e-6 else f"{conc*1e9}nM",
        "concentration_M": conc,
        "folder_name": folder_name
    })

# 生成sample sheet
df = pd.DataFrame(data)
df.to_csv("config/samples/TF_samples.tsv", sep='\t', index=False)
print("Sample sheet created successfully!")
print(df)
```

## 向后兼容性

新版本仍然支持旧的配置方式（直接列出sample_files和concentrations）。您可以：

1. **立即迁移**：使用sample sheet，获得更好的配置体验
2. **逐步迁移**：保持旧配置，在添加新TF时使用sample sheet
3. **不迁移**：继续使用旧配置方式

workflow会自动检测并使用正确的配置方式。

## 故障排除

### 错误：Sample sheet文件未找到

```
Error: Sample sheet for TF not found: config/samples/TF_samples.tsv
```

**解决**：检查文件路径是否正确，相对于workflow根目录。

### 错误：Sample file not found

```
Error: Sample file not found: /path/to/sample.txt
```

**解决**：
1. 检查 `sample_base_path` 是否正确
2. 检查 `sample_file_pattern` 格式是否匹配实际文件结构
3. 检查sample sheet中的 `folder_name` 是否正确

### 错误：Missing column

```
Error: Sample sheet missing column 'concentration_M'
```

**解决**：确保sample sheet包含所有必需列：`sample_id`, `concentration_M`, `folder_name`

## 获取帮助

如有问题，请：
1. 查看 `config/samples/README.md`
2. 运行配置验证：`python workflow/scripts/utils/validate_config.py`
3. 查看主README.md中的配置说明

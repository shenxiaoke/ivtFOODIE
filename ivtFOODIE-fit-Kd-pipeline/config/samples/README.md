# Sample Sheet 目录

本目录用于存放sample sheet文件（TSV格式）。

## Sample Sheet 格式说明

Sample sheet是制表符分隔的TSV文件，包含以下必需列：

| 列名 | 描述 | 示例 |
|------|------|------|
| `sample_id` | 样本ID（用于生成结果文件名） | `ZhangXueFei-1` |
| `concentration` | 浓度（可读标注，仅用于记录） | `2uM` 或 `500nM` |
| `concentration_M` | 浓度的摩尔值（用于Kd拟合） | `2.0e-6` |
| `folder_name` | 样本文件夹名称 | `251204-cA-mesc-T1-2uM_L3_NUDI001` |

## 创建Sample Sheet

### 方法1：从模板创建

复制 `template.tsv` 并填写您的样本信息：

```bash
cp template.tsv YourTF_samples.tsv
# 编辑 YourTF_samples.tsv
```

### 方法2：使用Excel/LibreOffice

1. 在Excel中创建表格，包含上述列
2. 另存为"制表符分隔的文本文件"（.txt或.tsv）
3. 确保编码为UTF-8

### 方法3：使用Python脚本生成

```python
import pandas as pd

samples = [
    {"sample_id": "TF-1", "concentration": "2uM", "concentration_M": 2.0e-6,
     "folder_name": "sample1_2uM"},
    {"sample_id": "TF-2", "concentration": "1uM", "concentration_M": 1.0e-6,
     "folder_name": "sample2_1uM"},
    # ... 添加更多样本
]

df = pd.DataFrame(samples)
df.to_csv("MyTF_samples.tsv", sep='\t', index=False)
```

## 在config.yaml中使用

创建sample sheet后，在 `config/config.yaml` 中配置：

```yaml
tfs:
  YourTF:
    species: "mouse"
    cell_line: "mESC"
    pc_control: "/path/to/PC_control_sites.txt"

    # Sample sheet配置
    sample_sheet: "config/samples/YourTF_samples.tsv"
    sample_base_path: "/path/to/samples/base/directory"
    sample_file_pattern: "{folder_name}/{folder_name}_sites.txt"

    fitting_indices: [0, 1, 2, 3, 4, 5]
    fix_n_to_1: true
```

`pc_control` 和 `sample_file_pattern` 都可以使用 `.txt` 或 `.txt.gz` 文件。`.txt.gz` 必须已经由bgzip压缩，并在同目录下带有对应的 `.txt.gz.tbi` 索引。例如：

```yaml
    pc_control: "/path/to/PC_control_sites.txt.gz"
    sample_file_pattern: "{folder_name}/{folder_name}_sites.txt.gz"
```

## 示例文件

- `ZhangXueFei_samples.tsv`: ZhangXueFei TF的样本配置
- `template.tsv`: 空白模板文件

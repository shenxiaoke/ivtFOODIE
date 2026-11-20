# 指定 Config 文件运行 Snakemake

## 基本用法

使用 `--configfile` 参数指定配置文件：

```bash
snakemake --configfile config/config-aid.yaml --cores 8
```

## 可用的配置文件

- `config/config-aid.yaml`
- `config/config-ssb.yaml`
- `config/config-aid_fix-ssb.yaml`
- `config/config-ssb_fix-aid.yaml`

## 常用命令

```bash
# 干运行（验证配置）
snakemake --configfile config/config-aid.yaml -n

# 完整运行
snakemake --configfile config/config-aid.yaml --cores 8

# 只运行到特定步骤
snakemake --configfile config/config-aid.yaml all_kd --cores 8

# 使用绝对路径
snakemake --configfile /absolute/path/to/config-aid.yaml --cores 8
```

## 注意事项

- 配置文件路径相对于运行 snakemake 命令的当前工作目录
- 建议在 `workflow/` 目录下运行，或使用相对路径 `../config/config-aid.yaml`
- 结果保存路径由 config 文件中的 `output.output_dir` 配置项指定

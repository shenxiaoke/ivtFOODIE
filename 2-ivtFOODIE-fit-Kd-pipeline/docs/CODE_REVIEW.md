# Snakemake 重构代码检查报告

## 阻塞性问题
- `workflow/rules/integrate.smk:12-54` 调用 `workflow/scripts/integration/integrate_samples.py` 作为 Snakemake `script`，但该脚本仅在 `__main__` 中通过 argparse 读取命令行参数（`...integrate_samples.py:379-415`），没有读取 `snakemake` 对象。Snakemake 调用时不会提供这些参数，会直接触发 argparse 的必填参数错误并导致步骤失败。
- 同一脚本的样本发现逻辑与当前工作目录结构不兼容。`integrate_samples_from_config` 假设 `sample_path` 下有按样本分组的子目录并在其中再查找 `.gz` 文件（`integrate_samples.py:317-347`），但预处理规则将压缩产物直接放在 `results/{tf}/01_compressed/*.gz`（`integrate.smk:45-54`）。对文件调用 `os.listdir(sample_folder)` 会抛出 `NotADirectoryError`，整合步骤无法运行。
- Flank 校正依赖的 motif 文件格式不匹配。`KdFitter._compute_flank_conversion` 强制读取 BED 的第 5 列作为 score（`kd_fitter.py:182-190`），默认回退路径 `results/{tf}/04_footprints/{tf}.footprint.basic.{thres}.bed` 来自 `denovo_footprint.py`，该输出只有 3 列（无 score）。在未提供打分 motif 时会产生 IndexError，flank 校正 Kd 拟合无法执行。

## 其他风险与缺口
- 迁移计划中标记的 motif 工具与可视化（`MIGRATION.md` 待实现项）尚未落地，但 README 与配置文件已暴露相关参数（如 STREME/FIMO、motif 可视化），用户可能误以为可用。
- 核心路径（如参考基因组、开放区域 BED）在 `config/config.yaml` 中硬编码为集群绝对路径，缺乏存在性校验或更友好的错误提示，易在新环境直接失败。

## 建议
1. 将 `integrate_samples.py` 改为读取 `snakemake.input/params/output/log`，或在 `rule integrate_samples` 中改用 `shell`/`script:` 传参方式，确保 DAG 能运行。顺便去掉对目录扫描的依赖，直接消费 Snakemake 已解析好的输入文件列表。
2. 为 flank 校正接受无 score 的 BED：允许 3 列输入时跳过打分筛选或使用均一得分；或在 `Snakefile` 中要求显式提供含 score 的 motif 文件，避免默认退回 footprint 结果。
3. 在配置加载或规则前置检查中添加文件存在性验证与更清晰的报错，减少运行期才发现路径无效的情况。

# 潜在缺陷总结

面向当前 `ivtFOODIE_postprocess` workflow 的一次静态审查，未在本地执行任何 `Snakemake` 任务或依赖计算。以下问题均可在现有代码上直接复现，请在远端正式环境中修复与验证。

## 1. Footprint 合并逻辑丢失重叠区间
- **位置**：`workflow/scripts/footprinting/denovo_footprint.py:247-263`
- **问题**：`process_rows_parallel()` 中声称“合并重叠 footprints”，但实现只是检测到重叠时 `continue`，直接丢弃后续片段，保留较早的 footprint。对于转录因子（Transcription Factor，TF）在高密度结合区形成的多个高 `Kd`（dissociation constant，解离常数）窗口，这会把真实事件截断为单个短片段，后续 motif 发现与可视化都将漏掉边界信息。
- **建议**：真正执行合并逻辑，例如在检测到重叠时更新 `prev_ed = max(prev_ed, ed)`，并按 `min(prev_st, st)`/`max(prev_ed, ed)` 写回，或直接调用 `bedtools merge` 保证区间完整。

## 2. Flank 校正无法处理空或低质量 motif
- **位置**：`workflow/scripts/fitting/kd_fitter.py:183-279`
- **问题**：
  1. 当 preliminary motif（默认来源于 basic footprint BED（Browser Extensible Data，浏览器可扩展数据））文件为空时，`motifs[0]` 会触发 `IndexError`；若只有极少 number 的 motif，`np.random.choice` 在 `len(motifs)=0` 时也会抛异常。
  2. 即便跳过异常，在没有有效 flank profile 时，`np.nanmean()` 返回 `NaN`，导致 `flank_baselines` 为 `NaN`，进一步让所有 `ratio_i_flank_adjust` 变成 `NaN` 并被 `dropna()` 整体过滤，`fit_kd_flank` 会输出空表。
- **影响**：对于信号弱或 footprint 为零的 TF，整个 flank corrected 分支直接崩溃。
- **建议**：在读入 motif 之后立即检查 `len(motifs)`，若为 0 则给出可理解的错误信息并回退到 basic 拟合；在计算 `flank_baselines` 时设置最小 motif 数量阈值，并对 `NaN` baseline 回落到 1.0 或跳过 flank 校正。

## 3. 固定 $n=1$ 时缺乏 Kd 取值约束
- **位置**：`workflow/scripts/fitting/kd_fitter.py:74-110`
- **问题**：`fix_n_to_1=True` 分支调用 `scipy.optimize.curve_fit` 时没有设置 `bounds`。一旦噪声较大，拟合可能产生负的 `Kd`，随后：
  - `hill_equation_n1()` 计算出的结合分数会超出 $[0,1]$。
  - `deltaG = -RT \ln(Kd)` 会因为 `Kd < 0` 抛出 `RuntimeWarning` 或返回复数。
  - `denovo_footprint` 阶段对 $-\log_{10}(Kd)$ 的分位数统计将出现 `NaN/inf`，直接破坏阈值计算。
- **建议**：与 `fix_n_to_1=False` 分支一致，为 `curve_fit` 添加 `bounds=(0, np.inf)`，或在拟合结果后强制 `kd = max(kd, 1e-15)`，保证物理意义成立。

## 4. `pc_name` 配置项未真正生效
- **位置**：`workflow/scripts/integration/integrate_samples.py:253-262` 以及 `workflow/rules/fitting.smk:32-43`
- **问题**：配置文件允许修改 `pc_name`（Protein Control，对照样本列名），`fit_kd_basic/fit_kd_flank` 也把该参数传给 `kd_fitter`。然而 `RatioIntegrator._dict_to_dataframe()` 始终将列命名为 `ratio_adjust/depth_adjust`。只要用户在 `config/config.yaml` 中把 `pc_name` 改成其他名字，`kd_fitter` 就会因为找不到对应列直接抛 `KeyError`。
- **建议**：在 integration 阶段即读取 `config["pc_name"]`，输出列名与配置保持一致；或者在 `Snakefile` 中强制 `pc_name` 固定并移除该配置项，避免产生“可配置但无效”的假象。

---

> ⚠️ **验证提示**：上述问题都会在 `snakemake -n` 层面暴露，但真正触发行路径需要远端实际数据（例如空 footprint、修改 `pc_name` 等）。请在目标集群上按需构造最小化示例，以免浪费正式批处理资源。

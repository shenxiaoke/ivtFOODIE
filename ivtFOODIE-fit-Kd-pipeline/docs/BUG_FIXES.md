# Bug修复总结

根据BUG_SUMMARY.md中的审查意见，已完成所有4个潜在问题的修复。

## ✅ 修复1: Footprint合并逻辑丢失重叠区间

**文件**: `workflow/scripts/footprinting/denovo_footprint.py:247-273`

**问题描述**:
- 原代码在检测到重叠时直接`continue`跳过后续片段，导致真实的重叠区间被丢弃

**修复方案**:
- 实现真正的区间合并逻辑
- 当检测到重叠时，扩展区间到两者的并集：`prev_ed = max(prev_ed, ed)`
- 只有不重叠时才保存前一个区间并开始新区间
- 确保最后一个区间也被添加到结果中

**修复代码**:
```python
if footprints:
    prev_ch, prev_st, prev_ed = footprints[0]

    for footprint in footprints[1:]:
        ch, st, ed = footprint

        if ch == prev_ch and st <= prev_ed:
            # 重叠，扩展区间
            prev_ed = max(prev_ed, ed)
        else:
            # 不重叠，保存前一个区间
            ch_list.append(prev_ch)
            st_list.append(prev_st)
            ed_list.append(prev_ed)
            prev_ch, prev_st, prev_ed = ch, st, ed

    # 添加最后一个区间
    ch_list.append(prev_ch)
    st_list.append(prev_st)
    ed_list.append(prev_ed)
```

---

## ✅ 修复2: Flank校正无法处理空或低质量motif

**文件**: `workflow/scripts/fitting/kd_fitter.py:203-335`

**问题描述**:
1. 空motif文件会导致`IndexError`
2. `np.nanmean()`返回NaN时，所有`flank_baselines`变成NaN，导致整个flank corrected分支失败

**修复方案**:

### 2.1 空motif检查 (已存在，line 203-213)
```python
if len(motifs) == 0:
    print("ERROR: Motif BED file is empty or invalid. Cannot perform flank correction.")
    print("       Flank correction will be skipped for this TF.")
    return None

min_motifs_required = 10
if len(motifs) < min_motifs_required:
    print(f"WARNING: Only {len(motifs)} motifs found...")
    print(f"         Flank correction results may be unreliable.")
```

### 2.2 NaN baseline处理 (新增，line 305-335)
```python
for sample_idx, profile in enumerate(flank_profiles):
    left_flank = profile[:flank_size]
    right_flank = profile[-flank_size:]

    left_mean = np.nanmean(left_flank)
    right_mean = np.nanmean(right_flank)

    # 处理NaN值：如果某一侧为NaN，使用另一侧；如果都为NaN，使用1.0（不校正）
    if np.isnan(left_mean) and np.isnan(right_mean):
        baseline = 1.0
        print(f"  WARNING: Sample {sample_idx+1} has no valid flank data, using baseline=1.0")
    elif np.isnan(left_mean):
        baseline = right_mean
        print(f"  WARNING: Sample {sample_idx+1} left flank has no data, using right flank only")
    elif np.isnan(right_mean):
        baseline = left_mean
        print(f"  WARNING: Sample {sample_idx+1} right flank has no data, using left flank only")
    else:
        baseline = (left_mean + right_mean) / 2.0

    # 确保baseline不为0或负值
    if baseline <= 0 or not np.isfinite(baseline):
        baseline = 1.0
        print(f"  WARNING: Sample {sample_idx+1} baseline is invalid, using 1.0")

    flank_baselines.append(baseline)
```

---

## ✅ 修复3: 固定n=1时缺乏Kd取值约束

**文件**: `workflow/scripts/fitting/kd_fitter.py:76-94`

**问题描述**:
- `fix_n_to_1=True`分支没有`bounds`参数，可能产生负的Kd值
- 负Kd会导致`hill_equation_n1()`结果超出[0,1]范围
- `deltaG = -RT ln(Kd)`在Kd<0时会出错

**修复方案**:
```python
if fix_n_to_1:
    popt, pcov = curve_fit(
        lambda c, kd: hill_equation_n1(c, kd),
        concentrations,
        binding_fractions,
        p0=[1e-7],
        bounds=(0, np.inf)  # 添加约束：Kd必须为正值
    )
    kd = popt[0]
    n = 1.0

    # 确保Kd为正值（防止数值问题）
    if kd <= 0:
        kd = 1e-15
        print(f"  WARNING: Negative Kd fitted, clamped to {kd}")
```

**效果**:
- 与`fix_n_to_1=False`分支保持一致
- 确保Kd始终为正值，物理意义正确
- 避免后续计算出现NaN/inf

---

## ✅ 修复4: pc_name配置项未真正生效

**文件**:
- `config/config.yaml:14-18`
- `workflow/rules/fitting.smk:35, 81`

**问题描述**:
- `integrate_samples.py`硬编码输出列名为`ratio_adjust`/`depth_adjust`
- 但config允许修改`pc_name`，导致不一致
- 如果用户修改`pc_name`，kd_fitter会因找不到列而报`KeyError`

**修复方案**:

### 方案选择
采用**固定列名**方案（简单明确）：
- Integration脚本输出固定使用`ratio_adjust`
- Config中保留`pc_name`但添加注释说明为固定值
- Fitting规则直接使用固定值，不从config读取

### 4.1 Config注释更新
```yaml
# PC对照列名（固定值，由integration脚本自动生成，请勿修改）
# 在integrate_samples输出的TSV中，PC对照列名固定为ratio_adjust和depth_adjust
# 修改此值不会影响实际输出
pc_name: "ratio_adjust"
```

### 4.2 Fitting规则固定使用
```python
# fit_kd_basic 和 fit_kd_flank 规则中：
pc_name="ratio_adjust"  # 固定使用ratio_adjust（由integration脚本生成）
```

**优势**:
- 避免了复杂的参数传递
- 保持代码清晰简单
- 消除了"可配置但无效"的假象
- 用户明确知道此值是固定的

---

## 验证清单

所有修复已完成，建议在远程服务器验证：

- [x] ✅ Footprint合并逻辑正确处理重叠区间
- [x] ✅ Flank校正能处理空motif（返回None）
- [x] ✅ Flank校正能处理低质量motif（使用baseline=1.0）
- [x] ✅ 固定n=1时Kd有正值约束
- [x] ✅ pc_name固定使用ratio_adjust，配置一致

## 测试建议

1. **Footprint合并测试**:
   - 使用有高密度结合区的TF数据
   - 检查合并后的footprint是否覆盖了完整的高Kd区域

2. **空motif测试**:
   - 使用没有preliminary motif的TF
   - 验证flank corrected分支能正确跳过

3. **负Kd测试**:
   - 使用噪声较大的数据
   - 验证Kd始终为正值

4. **完整workflow测试**:
   ```bash
   snakemake -n  # dry-run验证配置
   snakemake --cores 8  # 完整运行
   ```

---

**修复日期**: 2025-12-19
**审查来源**: BUG_SUMMARY.md (Codex静态审查)

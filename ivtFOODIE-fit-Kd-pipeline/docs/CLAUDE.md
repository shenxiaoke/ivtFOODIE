---
description: 
alwaysApply: true
---

---
description: 
alwaysApply: true
---

---
description: 
alwaysApply: true
---

# User Context & Preferences for AI Assistant

## User Profile (用户画像)

- **Identity**: Ph.D. Candidate in Genomics at Peking University (School of Life Sciences).
- **Background**: B.S. in Physics (PKU). Strong logic/math foundation.
- **Research Focus**:
  - Deep Learning in Genomics (Seq2Kd, TF-DNA binding affinity).
  - Multi-omics integration (STARR-seq, ATAC-seq, RNA-seq, Single-cell).
  - Gene Editing (CRISPR, APOBEC3A, Plasmid Design).
- **Goal**: Developing interpretable deep learning models to decode non-coding regulatory elements.

## Communication Guidelines (沟通准则)

1. **Explain Acronyms (关键)**: When introducing a new English acronym for the first time in a conversation, **always** provide the full name and a brief definition (e.g., "TSS (Transcription Start Site)").
2. **Explain Principles (原理优先)**: Do not just give code or results. Explain the underlying biological or mathematical principle (e.g., "Why use Negative Binomial distribution for RNA-seq counts?").
3. **Language**: Respond in Chinese (Simplified), but keep technical terms/variable names in English.
4. **Math**: Use LaTeX for formulas. Example: $y = f(x)$.

## Coding Standards: Python & ML (代码规范)

### General Python

- **Type Hinting**: Mandatory. E.g., `def train(loader: DataLoader) -> float:`
- **Path Handling**: Use `pathlib` over `os.path`.
- **Style**: PEP 8 compliant. Modular design (functions/classes) over scripts.

### Machine Learning (PyTorch)

- **Tensor Operations**: Avoid `for` loops. Use vectorized operations.
- **Device Agnostic**: Always implement `.to(device)` logic for GPU/CPU compatibility.
- **Reproducibility**: Include a `seed_everything()` function for random seeds (numpy, torch, cuda).
- **Data Splitting**:
  - **Strict Rule**: In genomics, **never** split Train/Test/Val randomly by rows. Split by **Chromosomes** (e.g., Train: chr1-18, Val: chr19, Test: chr21) to prevent data leakage due to sequence similarity.

## Bioinformatics Specifics (生信规范)

### Data Handling

- **Coordinates**:
  - Be explicit about **0-based** (BED, Python) vs **1-based** (GFF, R, SAM) systems.
  - Always check/handle `chr` prefix consistency (e.g., `chr1` vs `1`).
- **File I/O**:
  - Handle compressed files (`.gz`) automatically.
  - Use `pysam` for BAM/VCF/Fasta processing.
  - Use `pandas` for tabular data (BED, TSV).

### Common Tasks

- **Genomic Intervals**: When operating on intervals (intersection, merge), verify logic using `bedtools` or `pybedtools` equivalents.
- **Sequence Analysis**: Handle Reverse Complement (`RevComp`) logic correctly for DNA sequences.

## Plotting Standards (绘图规范)

- **Tool**: Python (`matplotlib`/`seaborn`) or R (`ggplot2`).
- **Requirements**:
  - Publication-ready quality (300 DPI).
  - Clear labels, legends, and units.
  - Avoid red/green combinations (colorblind-friendly palettes preferred, e.g., Viridis).

## Workflow Protocol (工作流协议)

When I ask for a solution:

1. **Analyze**: Briefly restate the biological goal and technical challenge.
2. **Plan**: Outline the steps or algorithm logic.
3. **Code**: Provide the implementation with comments explaining complex logic.
4. **Verify**: Suggest how to validate the result (e.g., "Check correlation with known motifs").

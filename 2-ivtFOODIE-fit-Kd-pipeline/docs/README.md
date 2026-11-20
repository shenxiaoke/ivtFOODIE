# ivtFOODIE postprocessing workflow

[中文版](README_CN.md)

This Snakemake workflow takes the per-site conversion ratio files from the mapping pipeline and fits dissociation constants (Kd) for individual cytosines. The full workflow also calls de novo footprints, discovers motifs with STREME, scans the genome with FIMO, and plots conversion ratios across concentrations at motif sites.

## Inputs

- Per-site conversion ratio files for each transcription factor (TF) at each tested protein concentration, with concentrations in molar units.
- A no-protein control (PC) ratio file for each TF. Set its path in `pc_control`.
- A reference genome FASTA for each species and an open-chromatin BED file for each cell line.
- An optional preliminary motif file for methods that use a known motif.

Provide ratio files as plain `.txt` or BGZF-compressed `.txt.gz`. The workflow compresses and indexes plain text inputs with `bgzip` and `tabix`. Precompressed inputs require a matching `.txt.gz.tbi` index in the same directory.

## Configure the samples

Start with a [configuration example](../config/config-SP1.yaml). Set the reference, open regions, sample and PC paths, and `output.output_dir` for your environment.

Use a [sample sheet](../config/samples/README.md) when files share a naming pattern. Each row supplies a `sample_id` for output filenames, a numeric `concentration_M` for fitting, and a `folder_name` for constructing the input path. The optional `concentration` column is a human-readable label. For example:

```yaml
tfs:
  YourTF:
    species: mouse
    cell_line: mESC
    pc_control: /path/to/PC_control_sites.txt
    sample_sheet: config/samples/YourTF.tsv
    sample_base_path: /path/to/ratio_files
    sample_file_pattern: "{folder_name}/{folder_name}_sites.txt"
    fix_n_to_1: true
```

When files do not share a naming pattern, specify `sample_files` and `concentrations` directly under the TF, with both lists in the same order. If you set `fitting_indices`, use zero-based positions in the sample sheet or file list.

## Conversion quality control

By default, the workflow compresses and indexes the ratio files, then plots mean conversion versus concentration before sample integration and Kd fitting. The PC is plotted at zero concentration. To compare multiple TFs in the same plot, include them in the same YAML configuration.

```yaml
conversion_qc:
  enabled: true
  n_sites: 10000
  sample_regions: 3000
  min_depth: 1
  seed: 261111
  mem_mb: 16000
```

For each TF, the workflow samples sites within randomly selected, merged open regions. A site must meet the QC depth threshold in every concentration sample and in the PC. All points on a TF curve use the same sites. Each plotted value is the unweighted arithmetic mean of the per-site conversion ratios, without PC normalization. If there are fewer eligible sites than requested, the workflow uses all available sites and reports a warning. QC fails if there are no eligible sites.

The QC depth threshold is independent of the top-level `min_depth` used for sample integration. Set `conversion_qc.enabled: false` to skip QC.

QC writes its results to `00_qc/<TF-group>/` under the configured output directory. The results include a PNG, an editable PDF, `mean_conversion.tsv` with per-sample mean and median depths and the fraction of sites with depth at least 10, a sampled-sites TSV for each TF, and `sampling_metadata.json` with sampling settings and the source PC paths.

## Run the workflow

From the Kd fitting pipeline directory, validate the selected configuration and preview the full workflow.

```bash
python workflow/scripts/utils/validate_config.py --config /absolute/path/config.yaml
snakemake --snakefile workflow/Snakefile --configfile /absolute/path/config.yaml --dry-run
```

Run locally with the workflow dependencies available in your environment.

```bash
snakemake --snakefile workflow/Snakefile --configfile /absolute/path/config.yaml --cores 8
```

To run only through a particular stage, append one of these targets to the Snakemake command.

| Target | Scope |
| --- | --- |
| `all_compressed` | Sample compression and indexing |
| `all_qc` | Conversion QC and its required input preparation |
| `all_integrated` | Sample integration |
| `all_kd` | Kd fitting |

For Slurm, use [`submit_snakemake.sh`](../submit_snakemake.sh) for a single configuration or [`submit_all_configs.sh`](../submit_all_configs.sh) for multiple configurations. The submission wrapper expects the xielab `bismark` Conda environment. Review the environment path and the configuration's cluster settings before using it on another cluster.

## Full workflow outputs

- Fitted Kd values for individual cytosine sites.
- De novo footprint sets called at the configured Kd thresholds.
- STREME motifs enriched in each footprint set.
- FIMO genomic matches for the enriched motifs.
- Concentration-dependent conversion ratios at motif sites.

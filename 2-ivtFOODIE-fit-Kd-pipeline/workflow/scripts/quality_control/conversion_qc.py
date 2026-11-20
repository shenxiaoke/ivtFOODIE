"""Reproducible, matched-site conversion/depth QC before sample integration."""
import argparse
from contextlib import ExitStack
from pathlib import Path
import json
import warnings

import numpy as np
import pandas as pd
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter


def merged_regions(path):
    intervals = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            chrom, start, end, *_ = line.split()
            start, end = int(start), int(end)
            if start < 0 or end <= start:
                raise ValueError('Invalid BED interval: ' + line)
            intervals.append((chrom, start, end))
    result = []
    for chrom, start, end in sorted(intervals):
        if result and result[-1][0] == chrom and start <= result[-1][2]:
            result[-1] = (chrom, result[-1][1], max(end, result[-1][2]))
        else:
            result.append((chrom, start, end))
    if not result:
        raise ValueError('Open-region BED is empty: ' + str(path))
    return result


def read_region(handle, chrom, start, end, min_depth):
    if chrom not in handle.contigs:
        return {}
    result = {}
    # Existing workflow indexes point coordinates without -0. Fetch one extra
    # base on either side and filter explicit 0-based positions afterwards.
    for line in handle.fetch(chrom, max(0, start - 1), end + 1):
        fields = line.split('\t')
        pos = int(fields[1])
        if not start <= pos < end:
            continue
        converted, unconverted = int(fields[3]), int(fields[4])
        depth = converted + unconverted
        if converted < 0 or unconverted < 0:
            raise ValueError('Negative conversion counts')
        if depth < min_depth:
            continue
        ratio = float(fields[5])
        if not np.isfinite(ratio) or not 0 <= ratio <= 1 or not np.isclose(ratio, converted / depth, atol=1e-5):
            raise ValueError('Invalid conversion ratio at ' + chrom + ':' + str(pos))
        if pos in result:
            raise ValueError('Duplicate ratio position at ' + chrom + ':' + str(pos))
        result[pos] = (ratio, depth)
    return result


def run_qc(spec, output_dir):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    options = spec.get('options', {})
    n_sites = int(options.get('n_sites', 10000))
    min_depth = int(options.get('min_depth', 1))
    seed = int(options.get('seed', 261111))
    n_regions = int(options.get('sample_regions', 3000))
    if min(n_sites, min_depth, n_regions) < 1:
        raise ValueError('n_sites, min_depth and sample_regions must be positive')
    rng = np.random.default_rng(seed)
    region_cache, summaries, audits = {}, [], []
    with ExitStack() as stack:
        opened = {}
        def open_tabix(path):
            path = str(Path(path).resolve())
            if path not in opened:
                opened[path] = stack.enter_context(pysam.TabixFile(path))
            return opened[path]

        for tf, info in sorted(spec['tfs'].items()):
            bed = info['open_regions']
            if bed not in region_cache:
                all_regions = merged_regions(bed)
                selected = rng.choice(len(all_regions), min(n_regions, len(all_regions)), replace=False)
                regions = sorted(all_regions[i] for i in selected)
                region_cache[bed] = regions
                with (out / ('sampled_regions_' + str(len(region_cache)) + '.bed')).open('w') as handle:
                    for chrom, start, end in regions:
                        handle.write(f'{chrom}\t{start}\t{end}\n')
            regions = region_cache[bed]
            entries = [dict(sample_id='PC', concentration_M=0., path=info['pc_control'], original_path=info.get('pc_source',info['pc_control']))] + info['samples']
            ids = [entry['sample_id'] for entry in entries]
            if len(ids) != len(set(ids)):
                raise ValueError('Sample IDs must be unique and must not equal PC: ' + tf)
            handles = [open_tabix(entry['path']) for entry in entries]
            keys, blocks = [], []
            for chrom, start, end in regions:
                data = [read_region(handle, chrom, start, end, min_depth) for handle in handles]
                common = set(data[0])
                for values in data[1:]:
                    common.intersection_update(values)
                positions = sorted(common)
                if positions:
                    keys.extend((chrom, pos) for pos in positions)
                    blocks.append(np.array([[values[pos] for values in data] for pos in positions], dtype=float))
            if not keys:
                raise ValueError('No shared depth-qualified sites in sampled regions for ' + tf)
            count = min(n_sites, len(keys))
            if count < n_sites:
                warnings.warn(f'{tf}: only {count} common sites; requested {n_sites}. See metadata.')
            indices = rng.choice(len(keys), count, replace=False)
            chosen = [keys[i] for i in indices]
            values = np.concatenate(blocks)[indices]
            site_table = pd.DataFrame(dict(chrom=[k[0] for k in chosen], start=[k[1] for k in chosen]))
            site_table['end'] = site_table['start'] + 1
            for j, entry in enumerate(entries):
                ratios, depths = values[:, j, 0], values[:, j, 1]
                site_table[entry['sample_id'] + '_ratio'] = ratios
                site_table[entry['sample_id'] + '_depth'] = depths.astype(int)
                summaries.append(dict(tf=tf, sample_id=entry['sample_id'], source_path=entry.get('original_path',entry['path']), concentration_uM=float(entry['concentration_M']) * 1e6, n_sites=count, mean_conversion=float(ratios.mean()), mean_conversion_percent=float(ratios.mean()*100), mean_depth=float(depths.mean()), median_depth=float(np.median(depths)), fraction_depth_ge10=float((depths>=10).mean())))
            site_table.to_csv(out / (tf + '.sampled_sites.tsv'), sep='\t', index=False)
            audits.append(dict(tf=tf, candidate_common_sites=len(keys), sampled_sites=count, sampled_regions=len(regions), pc_control=info['pc_control']))
            print(f'{tf}: sampled {count:,}/{len(keys):,} shared eligible sites', flush=True)
    frame = pd.DataFrame(summaries).sort_values(['tf', 'concentration_uM'])
    frame.to_csv(out/'mean_conversion.tsv', sep='\t', index=False)
    render(frame, out, min_depth)
    metadata = dict(seed=seed, n_sites_requested=n_sites, min_depth=min_depth, audits=audits, inputs=spec, method='Random merged open-region clusters, then uniform sampling without replacement of common sites per TF. Within each TF, all samples and PC use the same sites. Mean is the unweighted arithmetic mean of site conversion rates; no PC normalization. Depth cutoff applies to QC only, independently of fitting min_depth.')
    (out/'sampling_metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    return frame


def render(frame, output_dir, min_depth):
    colors = ['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9','#882255','#332288','#666666','#88A825']
    markers = ['o','s','^','D','v','P','X','<','>','h']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'ps.fonttype':42})
    fig, ax = plt.subplots(figsize=(11.5,7.2))
    for i,(tf,group) in enumerate(frame.groupby('tf', sort=True)):
        ax.plot(group.concentration_uM,group.mean_conversion,color=colors[i%len(colors)],marker=markers[i%len(markers)],linewidth=1.8,markersize=5.8,label=tf,alpha=.94)
    ax.set(xlabel='Protein concentration (µM)',ylabel='Mean conversion rate',title='ivtFOODIE: mean conversion across protein concentrations',ylim=(0,1))
    concentrations = sorted(frame.concentration_uM.unique())
    # Numeric positioning preserves real dose spacing; thin dense tick labels.
    ticks = [concentrations[0]]
    gap = max(concentrations[-1] / 18, 1e-12)
    for value in concentrations[1:]:
        if value - ticks[-1] >= gap or value == concentrations[-1]:
            ticks.append(value)
    ax.set_xticks(ticks)
    ax.set_xticklabels(['0\n(PC)' if v==0 else f'{v:g}' for v in ticks])
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y',color='#dddddd',linewidth=.65);ax.set_axisbelow(True)
    ax.legend(title='TF',loc='upper left',bbox_to_anchor=(1.015,1),frameon=False,fontsize=10)
    pcs=frame[frame.sample_id=='PC'].source_path.map(lambda p:Path(p).name).unique()
    pc_note='Shared PC: '+pcs[0] if len(pcs)==1 else 'PC sources are recorded in mean_conversion.tsv.'
    fig.text(.09,.025,f'{frame.n_sites.min():,}–{frame.n_sites.max():,} matched sites per TF; QC depth ≥{min_depth}; unweighted means.\n{pc_note}\nPC can differ slightly across TFs because sampled site sets differ.',fontsize=8,color='#555555')
    fig.subplots_adjust(left=.09,right=.80,top=.90,bottom=.17)
    for ext in ['png','pdf']:
        fig.savefig(Path(output_dir)/('mean_conversion_vs_concentration.'+ext),dpi=220)
    plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spec',required=True,help='JSON specification (TF samples and indexed ratio paths)')
    parser.add_argument('--output-dir',required=True)
    args=parser.parse_args()
    run_qc(json.loads(args.spec),args.output_dir)

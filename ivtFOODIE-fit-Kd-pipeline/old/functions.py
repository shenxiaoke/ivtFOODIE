import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats


def create_sample_links():
    from_dir=/dshare/xielab/analysisdata/TFC/shenk/analysis/WZ/BANC_FOODIE/SP1_20241124/data/20241223-SP1ti-highhighYY1/02-merge-outputs/ratio
    target_dir=/dshare/xielab/analysisdata/TFC/shenk/analysis/WZ/TF_wholeGenome_inVitro_binding/TF_wholeGenome_inVitro_binding_20240923_batchAdjust_new/data/ivtFOODIE/SP1_BANCseq/02-merge-outputs/ratio
    
    old_values=(2d4u 2u 1d6u 1d2u 800n 600n 400n 200n 100n 50n)
    new_values=(1 2 3 4 5 6 7 8 9 10)
    
    # # 去掉测序深度浅的样本
    # old_values=(2d4u 800n 600n 400n 100n 50n)
    # new_values=(1 2 3 4 5 6)
    
    for i in "${!old_values[@]}"; do
        old="${old_values[i]}"
        new="${new_values[i]}"
        
        echo "Processing: $old -> $new"
    
        fin=`ls $from_dir|grep $old|head -1`
        fout=SP1_BANCseq-$new
    
        # echo $from_dir/$fin/${fin}_sites.txt.gz $target_dir/$fout/${fout}_sites.txt.gz
    
        mkdir -p $target_dir/$fout
        ln -s $from_dir/$fin/${fin}_sites.txt.gz $target_dir/$fout/${fout}_sites.txt.gz
        ln -s $from_dir/$fin/${fin}_sites.txt.gz.tbi $target_dir/$fout/${fout}_sites.txt.gz.tbi
    
    done
    pass

def generate_single_C_ratio_file():
    TF = 'PAX7_Da829'
    samples_from = 1
    samples_to = 7
    open_region_f = open_region_GM
    
    sample_path = f'/dshare/xielab/analysisdata/TFC/shenk/analysis/WZ/TF_wholeGenome_inVitro_binding/\
    TF_wholeGenome_inVitro_binding_20240923_batchAdjust_new/data/ivtFOODIE/{TF}/02-merge-outputs/ratio'
    
    PC_path = '/dshare/xielab/analysisdata/TFC/shenk/analysis/WZ/TF_wholeGenome_inVitro_binding/\
    TF_wholeGenome_inVitro_binding_20240923_batchAdjust_new/data/ivtFOODIE/PC/Sample_241030-GM-cA1017-0d1x-PC_sites.txt.gz'
    
    !mkdir -p ../results/C_site_wise_Kd_deduplicate_Da829
    out_path = '../results/C_site_wise_Kd_deduplicate_Da829'
    
    !sbatch -p compute_new -t 7-0:00 -c 1 --mem 16g --wrap="python write_samples_ratio.py --TF {TF} --sample_path {sample_path} --samples_from {samples_from} --samples_to {samples_to} --PC_path {PC_path} --open_region_f {open_region_f} --out_path {out_path}"

    pass

def fitting_Kd():
    # old_values=(2uM 1d5uM 800nM 500nM 250nM)
    # new_values=(1 2 3 4 5)
    
    in_f='../results/C_site_wise_Kd_deduplicate_ZhangXueFei/ZhangXueFei.samples_ratio.csv'
    out_f='../results/C_site_wise_Kd_deduplicate_ZhangXueFei/ZhangXueFei.fittingKd.n_equal_1.dropna.csv'
    concentration='2e-6,1.5e-6,0.8e-6,0.5e-6,0.25e-6'
    fitting_used_index='0,1,2,3,4' # from 1 to 5
    fix_n_to_1='True'
    
    sbatch -p compute_plus -c 1 --mem 16g -t 7-0:00 \
    --wrap="/share/home/shenk/local/anaconda3/envs/bismark/bin/python Kd_noBatchAdjust.py --in_f $in_f --out_f $out_f --concentration $concentration --fitting_used_index $fitting_used_index --fix_n_to_1 $fix_n_to_1"

    pass

def fitting_Kd_flankAdjust():
    # old_values=(2uM 1d5uM 800nM 500nM 250nM)
    # new_values=(1 2 3 4 5)
    
    in_f='../results/C_site_wise_Kd_deduplicate_ZhangXueFei/ZhangXueFei.samples_ratio.csv'
    out_f='../results/C_site_wise_Kd_deduplicate_ZhangXueFei/ZhangXueFei.fittingKd.n_equal_1.dropna.flank_adjust.csv'
    concentration='2e-6,1.5e-6,0.8e-6,0.5e-6,0.25e-6'
    fitting_used_index='0,1,2,3,4' # from 1 to 5
    fix_n_to_1='True'
    motif_f='../results/Tox_motif/Sp1.open.bed'
    
    sbatch -p compute_plus -c 1 --mem 16g -t 7-0:00 \
    --wrap="/share/home/shenk/local/anaconda3/envs/bismark/bin/python Kd_flankRatioCorrected.py --in_f $in_f --out_f $out_f --motif_f $motif_f --concentration $concentration --fitting_used_index $fitting_used_index --fix_n_to_1 $fix_n_to_1"

    pass

def de_novo_footprinting_by_automatically_selected_Kd_threshold():
    import os
    import pandas as pd
    import numpy as np
    
    HM='../results/C_site_wise_Kd_deduplicate_ZhangXueFei'
    out_hm = '../results/C_site_wise_Kd_deduplicate_ZhangXueFei_footprint_automatically_selected/top5_top10_Kd_threshold'
    # out_hm = '../results/C_site_wise_Kd_deduplicate_ZhangXueFei_footprint_automatically_selected/top10_top20_Kd_threshold'
    # out_hm = '../results/C_site_wise_Kd_deduplicate_ZhangXueFei_footprint_automatically_selected/top15_top30_Kd_threshold'
    !mkdir -p $out_hm
    
    exist_files = [i.split('.footprint')[0] for i in os.listdir(out_hm)]
    
    num_processes = 16
    
    for f in os.listdir(HM):
        if 'fittingKd' not in f or 'K562' in f or 'dddB' in f or 'flank_adjust' in f:
            continue
        # print(f)
        in_f = f'{HM}/{f}'
        f_name = f.split('.fittingKd')[0]
        
        df = pd.read_csv(in_f,sep='\t', nrows=10000)
    
        Kd_log = -np.log10(df['Kd'])
    
        deltaG_thres_seed, deltaG_thres = round(np.quantile(Kd_log,0.95),2), round(np.quantile(Kd_log,0.9),2)
        # deltaG_thres_seed, deltaG_thres = round(np.quantile(Kd_log,0.9),2), round(np.quantile(Kd_log,0.8),2)
        # deltaG_thres_seed, deltaG_thres = round(np.quantile(Kd_log,0.85),2), round(np.quantile(Kd_log,0.7),2)
        out_f = f'{out_hm}/{f_name}.footprint.{deltaG_thres_seed}_{deltaG_thres}.bed'
    
        if f_name in exist_files:
            continue
        print(f_name)
        !sbatch -p compute_plus -t 7-0:00 -c 16 --mem 32g --wrap="echo $in_f $out_f $deltaG_thres_seed $deltaG_thres $num_processes|~/local/anaconda3/envs/bismark/bin/python footprinting.py"
    pass

from scipy.special import xlogy

def compl(seq):
    comp = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
    return ''.join([comp[i] for i in seq])[::-1]

def read_pfm(file):
    info_list = []
    pfm_list = []
    with open(file) as f:
        for l in f:
            if 'MOTIF' in l:
                # motif_num = int(l.strip().split('-')[-1])
                motif_num = 1
                pfm = []
                info = next(f).strip()
                width,nsites = int(info.split('w= ')[1].split(' ')[0]),int(info.split('nsites= ')[1].split(' ')[0])
                info_list.append([motif_num,width,nsites])
                l = next(f)
                while ('.' in l): # non-empty line
                    pfm.append(list(map(float,l.strip().split())))
                    l = next(f)
                pfm_list.append(np.transpose(np.array(pfm),(1,0)))
    return info_list,pfm_list

def pfm2pwmlogo(pfm):
    pfm = pfm/np.sum(pfm,axis=0).reshape((1,-1))
    R = np.array(2-(-xlogy(pfm,pfm) * np.log2(np.e) ).sum(axis=0))
    R = R.reshape((1,-1))
    pwm = np.array(pfm)*R
    return pwm

# real pwm instead of pwm logo
def pfm2pwm(pfm):
    pfm = pfm/np.sum(pfm,axis=0).reshape((1,-1))
    epsilon = 1e-4
    pwm = np.log2((pfm + epsilon)/0.25)
    return pwm

def logo_plot(pfm,ax=None):
    df_pfm = pd.DataFrame(np.transpose(pfm,(1,0)))
    # df_pfm = pd.DataFrame(pfm)
    df_pfm.columns = ['A','C','G','T']

    if ax:
        logo = logomaker.Logo(df_pfm,
                              ax=ax,
                              width=.8,
                              vpad=.05,)
    else:
        logo = logomaker.Logo(df_pfm,
                              width=.8,
                              vpad=.05,)

    # style using Logo methods
    logo.style_spines(spines=['left', 'right', 'top'], visible=False)

    # style using Axes methods
    logo.ax.set_xticks(range(len(df_pfm)))
    # logo.ax.set_xticklabels('%+d'%x for x in [-3, -2, -1, 1, 2, 3, 4, 5, 6])
    logo.ax.set_yticks([0, .5, 1])
    # logo.ax.axvline(2.5, color='k', linewidth=1, linestyle=':')
    logo.ax.set_ylabel('probability')
    # plt.show()

def motif_scan():
    sbatch -p compute_plus -c 1 --mem 16g -t 2-0:00 \
--wrap="~/meme/bin/fimo --thresh 0.001 --max-stored-scores 100000000 --verbosity 5  --oc outputs_ZhangXueFeiAID_DGYW_maxStoredScores10M_p0d001  ZhangXueFeiAID_DGYW.meme /dshare/xielab/analysisdata/TFC/shenk/important_data/chrom/genome/Bismark_mm10/mm10.genome.fa"

    # pwm = pfm2pwm(pfm_list[i])
    # pwm
    # pwmlogo = pfm2pwmlogo(pfm_list[i])
    # logo_plot(pwmlogo)
    # plt.show()
    
    pass

class motifCenteredPlot_siteWise:
    def __init__(self,gz_f,gz_control_f,motif_f,motif_pick_ratio=0.1,span=100,min_depth=10):
        self.gz = tabix.open(gz_f)
        self.gz_control = tabix.open(gz_control_f)
        self.motif_f = motif_f
        self.motif_pick_ratio = motif_pick_ratio
        self.span = span
        self.min_depth = min_depth
    def motifs_list(self,pick_top_or_bottom='top',motif_pick_ratio=None):
        if motif_pick_ratio!=None:
            self.motif_pick_ratio = motif_pick_ratio
        motifs = []
        motif_score_list = []
        count = 0
        with open(self.motif_f,'r') as f:
            for l in f:
                count += 1
                if count % 1000 == 0:
                    print(count,end='\r')
                l = l.strip().split('\t')
                motifs.append(l)
                motif_score_list.append(float(l[4]))
        if pick_top_or_bottom=='top':
            self.score_thres = np.percentile(motif_score_list,100*(1-self.motif_pick_ratio))
            self.chosen_motifs = list(filter(lambda x:float(x[4])>=self.score_thres, motifs))
        elif pick_top_or_bottom=='bottom':
            self.score_thres = np.percentile(motif_score_list,100*(self.motif_pick_ratio))
            self.chosen_motifs = list(filter(lambda x:float(x[4])<=self.score_thres, motifs))
        else:
            print('pick_top_or_bottom para error')
        print(f'chosen motifs number: {len(self.chosen_motifs)}')
        print(f'score_thres: {self.score_thres}')
        return self.chosen_motifs
    def conversion_profile(self,method='median'):
        self.conversion_list = [[] for i in range(self.span)]
        count = 0
        for motif in self.chosen_motifs:
            count += 1
            if count % 1000 == 0:
                print(count,end='\r')
            ch,st,ed,motif_name,score,strand,motif_seq = motif
            # ch,st,ed,motif_name,score,motif_id,strand,motif_seq,_,__ = motif
            center = (int(st)+int(ed))//2
            try:
                records = self.gz.query(ch,center-self.span//2,center+self.span//2)
                records_control = self.gz_control.query(ch,center-self.span//2,center+self.span//2)
            except:
                continue
            # print('hello')
            records = {record[1]:record for record in records}
            records_control = {record[1]:record for record in records_control}
            for pos in set(records.keys())&set(records_control.keys()):
                record,record_control = records[pos],records_control[pos]
                depth,depth_control = int(record[3])+int(record[4]),int(record_control[3])+int(record_control[4])
                ratio,ratio_control = float(record[5]),float(record_control[5])
                if depth<self.min_depth or depth_control<self.min_depth or ratio_control==0:
                    continue
                index = int(pos)-center+self.span//2 if strand=='+' else center-int(pos)+self.span//2
                if 0<=index<self.span:
                    self.conversion_list[index].append(ratio/ratio_control)
        self.conversion = np.array([np.median(i) if len(i)>0 else np.nan for i in self.conversion_list]) if method=='median' else \
                            np.array([np.mean(i) if len(i)>0 else np.nan for i in self.conversion_list])
        return self.conversion
    def plot(self,left_right_pos=None,add_conversion=None,figsize=None,normalize=False,normalize_flank=False,plt_save=None,ylim=None,label=None):
        y1 = add_conversion
        if left_right_pos==None:
            x = np.arange(-self.span//2,self.span//2)
            y = self.conversion
        else:
            left_pos,right_pos = left_right_pos
            x = np.arange(left_pos,right_pos)
            y = self.conversion[self.span//2+left_pos:self.span//2+right_pos]
            if add_conversion is not None:
                y1 = {i:y1[add_conversion][self.span//2+left_pos:self.span//2+right_pos] for i in y1}
        if normalize:
            factor = np.mean(y)
            if normalize_flank:
                factor = np.mean(list(y[:self.span//2-20])+list(y[self.span//2+20:]))
            y = np.array(y)/factor
            if add_conversion is not None:
                factor = {i:np.mean(y1[i]) for i in y1}
                if normalize_flank:
                    factor = {i:np.mean(list(y1[i][:self.span//2-20])+list(y1[i][self.span//2+20:])) for i in y1}
                y1 = {i:np.array(y1[i])/factor[i] for i in y1}
        # if figsize is not None:
        #     plt.subplots(figsize=figsize)
        # else:
        #     plt.subplots(figsize=(8,5))
        plt.plot(x,y,label='y') if label is None else plt.plot(x,y,label=label)
        if add_conversion is not None:
            for i in y1:
                plt.plot(x,y1[i],label=str(i))
        # plt.legend()
        if ylim:
            plt.ylim(ylim)
        if plt_save:
            plt.savefig(plt_save)
        # plt.show()
        return

def convert_single_C_Kd_to_igvtdf_format():
    pass


gz_control_f = '/dshare/xielab/analysisdata/TFC/shenk/rawdata/20251118_WZ_ivtFOODIE/02-merge-outputs/ratio_rename/251103-cAmesc-Da829B4-0d1x-PC/251103-cAmesc-Da829B4-0d1x-PC_sites.txt.gz'

sample_hm = f'/dshare/xielab/analysisdata/TFC/shenk/analysis/WZ/TF_wholeGenome_inVitro_binding/\
TF_wholeGenome_inVitro_binding_20240923_batchAdjust_new/data/ivtFOODIE/ZhangXueFei/02-merge-outputs/ratio'
    
TF = 'ZhangXueFei'
motif_f = f'../results/ZhangXueFei_motif/ZhangXueFeiIvtFOODIE_GCCCCGCCC.open.bed'

fig,ax = plt.subplots(figsize=(10,7))

for f in sorted(glob.glob(f'{sample_hm}/*')):
    f = f.split('/')[-1]
    print(f)

    gz_f = f'{sample_hm}/{f}/{f}_sites.txt.gz'
    
    a = motifCenteredPlot_siteWise(gz_f,gz_control_f,motif_f,min_depth=5)
    
    _ = a.motifs_list(pick_top_or_bottom='top',motif_pick_ratio=0.01)
    a.conversion_profile(method='median')
    # print(a.conversion)
    a.plot(normalize=True,normalize_flank=True,figsize=(10,7),label=f.split('-')[-1])
_ = plt.legend()
























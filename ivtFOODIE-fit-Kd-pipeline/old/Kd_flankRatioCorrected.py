import argparse
import numpy as np
import tabix,pysam
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit
import statsmodels.api as sm

from scipy.optimize import least_squares

base_list = ['A','C','G','T','N']
compl_dict = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
def complementary(context):
    return ''.join([compl_dict[i] for i in context])[::-1]

def ratio_split(ratio_df):
    ratio_dict = {}
    for ch in ratio_df['ch'].unique():
        ratio_dict[ch] = ratio_df[ratio_df['ch']==ch]
    return ratio_dict
def flank_conversion(ratio_df,motif_f,motif_pick_ratio=0.1,span=100,PC_name='ratio_adjust'):
    ratio_df_dict = ratio_split(ratio_df)
    motifs = []
    motif_score_list = []
    with open(motif_f,'r') as f:
        for l in f:
            l = l.strip().split('\t')
            motifs.append(l)
            motif_score_list.append(float(l[4]))
    score_thres = np.percentile(motif_score_list,100*(1-motif_pick_ratio))
    chosen_motifs = list(filter(lambda x:float(x[4])>=score_thres, motifs))

    # sample_list = [i for i in ratio_df.columns if ('ratio_' in i and 'ratio_adjust' not in i)]
    sample_list = [i for i in ratio_df.columns if ('ratio_' in i and PC_name not in i)]
    conv_list = [[[] for _ in range(2*span)] for __ in range(len(sample_list))]
    for number,motif in enumerate(chosen_motifs):
        if number % 100 == 0:
            print('motif number:',number, end='\r')
        ch,st,ed = motif[0],int(motif[1]),int(motif[2])
        center = (st+ed)//2
        if ch not in ratio_df_dict:
            continue
        selected_df = ratio_df_dict[ch]
        selected_df = selected_df[(selected_df['st']>=st-span)&(selected_df['st']<ed+span)]
        for i,row in selected_df.iterrows():
            C_pos = row['st']
            rela_pos = C_pos-center+span
            PC_ratio = row[PC_name]
            if PC_ratio==0:
                continue
            if 0<=rela_pos<2*span:
                for i in range(len(sample_list)):
                    sample = sample_list[i]
                    conversion = row[sample]
                    conv_list[i][rela_pos].append(conversion/PC_ratio)
    for i in range(len(conv_list)):
        conv_list[i] = [np.nanmedian(conv_list[i][j]) for j in range(2*span)]
    return conv_list

def flank_conversion_adjust(ratio_df, flank_conversion_list,flank=50,PC_name='ratio_adjust'):
    adjust_conv_list = [(np.nanmean(i[:flank])+np.nanmean(i[-flank:]))/2 for i in flank_conversion_list]
    # sample_list = [i for i in ratio_df.columns if ('ratio_' in i and 'ratio_adjust' not in i)]
    sample_list = [i for i in ratio_df.columns if ('ratio_' in i and PC_name not in i)]
    for sample,adjust_conv in zip(sample_list,adjust_conv_list):
        ratio_df[f'{sample}_flank_adjust'] = ratio_df[sample]/adjust_conv
    return ratio_df

base_list = ['A','C','G','T','N']
compl_dict = {'A':'T','T':'A','C':'G','G':'C','N':'N'}

def complementary(context):
    return ''.join([compl_dict[i] for i in context])[::-1]

def hill_fit(x,Kd,n):
    return 1/(1+(Kd/x)**n)

def hill_fit_n1(x,Kd):
    return 1/(1+(Kd/x))

def Kd_fit_siteWise(x, concentration, fitting_used_index = np.arange(0,13),batchAdjust=True,PC_name='ratio_adjust'):
    if x['counter']%1000==0:
        print(f"processed {x['counter']}")
    # print(x['counter'],end='\r')
    if batchAdjust:
        y = np.array([x[f'ratio_{i}_flank_adjust'] for i in fitting_used_index+1])
    else:
        y = np.array([x[f'ratio_{i}'] for i in fitting_used_index+1])
    
    used_concentration = np.array(concentration)[fitting_used_index]
    spearman_correlation = -stats.spearmanr(used_concentration, y)[0]

    PC_ratio = x[PC_name]
    if PC_ratio==0:
        return [np.nan,np.nan,spearman_correlation,np.nan,np.nan]
    y /= PC_ratio
    
    y = 1-y
    y = np.clip(y,0,1)
    try:
        popt, pcov = curve_fit(hill_fit, used_concentration, y, p0=[0.1e-6,1],bounds=((0,0), (1,10)))
    except RuntimeError:
        # return np.nan
        return [np.nan,np.nan,spearman_correlation,np.nan,np.nan]
    Kd,n = popt
    exp_y = np.array(hill_fit(used_concentration,Kd,n))
    err = 100*np.mean((exp_y-y)**2)
    KdErr = np.sqrt(np.diag(pcov))[0]/Kd
    return [Kd,n,spearman_correlation,err,KdErr]

def Kd_fit_n1_siteWise(x, concentration, fitting_used_index = np.arange(0,13),batchAdjust=True,PC_name='ratio_adjust'):
    if x['counter']%1000==0:
        print(f"processed {x['counter']}")
    # print(x['counter'],end='\r')
    if batchAdjust:
        y = np.array([x[f'ratio_{i}_flank_adjust'] for i in fitting_used_index+1])
    else:
        y = np.array([x[f'ratio_{i}'] for i in fitting_used_index+1])
    
    used_concentration = np.array(concentration)[fitting_used_index]
    spearman_correlation = -stats.spearmanr(used_concentration, y)[0]

    PC_ratio = x[PC_name]
    if PC_ratio==0:
        return [np.nan,np.nan,spearman_correlation,np.nan,np.nan]
    y /= PC_ratio

    y = 1-y
    y = np.clip(y,0,1)
    try:
        popt, pcov = curve_fit(hill_fit_n1, used_concentration, y, p0=[0.1e-6])
    except RuntimeError:
        # return np.nan
        return [np.nan,1,spearman_correlation,np.nan,np.nan]
    Kd = popt[0]
    exp_y = np.array(hill_fit_n1(used_concentration,Kd))
    err = 100*np.mean((exp_y-y)**2)
    KdErr = np.sqrt(np.diag(pcov))[0]/Kd
    return [Kd,1,spearman_correlation,err,KdErr]


def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description='fit Kd')

    # Add argument
    parser.add_argument('--in_f', type=str, help='input ratio file')
    parser.add_argument('--out_f', type=str, help='output Kd file')
    parser.add_argument('--concentration', type=str, help='concentration split by comma')
    parser.add_argument('--fitting_used_index', type=str, help='index split by comma')
    parser.add_argument('--fix_n_to_1', type=bool, help='whether or not fix hill coefficient to 1')
    parser.add_argument('--motif_f', type=str, help='motif file path')
    parser.add_argument('--PC_name', type=str, default='ratio_adjust', help='PC name')

    # Parse the arguments
    args = parser.parse_args()

    in_f = args.in_f
    out_f = args.out_f
    concentration = np.array( list(map(float,args.concentration.split(','))) )
    fitting_used_index = np.array( sorted(list(map(int,args.fitting_used_index.split(',')))) )
    fix_n_to_1 = args.fix_n_to_1
    motif_f = args.motif_f
    PC_name = args.PC_name

    # Print the argument
    print(in_f, out_f, concentration, fitting_used_index, fix_n_to_1,motif_f, PC_name)

    # fit
    ratio_df = pd.read_csv(in_f,sep='\t')
    print('ratio_df load finished')
    ratio_df['counter'] = ratio_df.index

    print(motif_f)
    flank_conversion_list = flank_conversion(ratio_df,motif_f,PC_name=PC_name)
    flank_conversion_list = [np.array(flank_conversion_list[i]) + 1e-6 for i in range(len(flank_conversion_list))] # avoid 0
    print('flank conversion done')

    
    ratio_df = flank_conversion_adjust(ratio_df, flank_conversion_list=flank_conversion_list,PC_name=PC_name)

    ratio_df = ratio_df.dropna().copy()

    print('flank adjust finished')
    
    if fix_n_to_1:
        ratio_df['Kd_n_spearmanr_err_KdErr'] = ratio_df.apply(lambda x:Kd_fit_n1_siteWise(x, concentration=concentration, fitting_used_index=fitting_used_index, PC_name=PC_name), axis=1)
    else:
        ratio_df['Kd_n_spearmanr_err_KdErr'] = ratio_df.apply(lambda x:Kd_fit_siteWise(x, concentration=concentration, fitting_used_index=fitting_used_index, PC_name=PC_name), axis=1)

    print('Kd fitting finished')
    
    ratio_df['Kd'] = ratio_df.apply(lambda x:x['Kd_n_spearmanr_err_KdErr'][0], axis=1)
    ratio_df['n'] = ratio_df.apply(lambda x:x['Kd_n_spearmanr_err_KdErr'][1], axis=1)
    ratio_df['spearmanr'] = ratio_df.apply(lambda x:x['Kd_n_spearmanr_err_KdErr'][2], axis=1)
    ratio_df['err'] = ratio_df.apply(lambda x:x['Kd_n_spearmanr_err_KdErr'][3], axis=1)
    ratio_df['KdErr'] = ratio_df.apply(lambda x:x['Kd_n_spearmanr_err_KdErr'][4], axis=1)
    ratio_df.drop('Kd_n_spearmanr_err_KdErr',inplace=True, axis=1)
    R=8.31446261815324 # J·K−1·mol−1
    T=300
    ratio_df['deltaG'] = -R*T*np.log(ratio_df['Kd'])/1e3 # kJ
    ratio_df.dropna().to_csv(out_f,sep='\t',index=False, na_rep='NULL')

    print('done')
    

if __name__ == "__main__":
    main()

import pandas as pd
import numpy as np
import sys,pysam

import os
import multiprocessing
from functools import partial

def ratio_split(ratio_df):
    ratio_dict = {}
    for ch in ratio_df['ch'].unique():
        ratio_dict[ch] = ratio_df[ratio_df['ch']==ch]
    return ratio_dict

def deNovoFootprinting(index,ratio_df,ratio_df_dict,deltaG_thres_seed=None,deltaG_thres=None,min_len=10,max_len=30,max_gap=5,tolerance=0.1, Kd_column='Kd_log'):
    count = 0
    ch,st,ed = [ ratio_df.iloc[index,][i] for i in ['ch','st','ed'] ]
    selected_ratio_df = ratio_df_dict[ch]
    selected_ratio_df = selected_ratio_df[ (selected_ratio_df['st']>=st-max_len)&(selected_ratio_df['ed']<st+max_len) ]
    st_list,deltaG_list = list(selected_ratio_df['st']),np.array(selected_ratio_df[Kd_column])
    center_i = st_list.index(st)
    stop_left,stop_right = False,False
    dormant_left,dormant_right = False,False
    left_i,right_i = center_i,center_i
    while stop_left==False or stop_right==False:
        if stop_left==False:
            if left_i-1<0:
                stop_left = True
            elif st_list[left_i]-st_list[left_i-1]>max_gap:
                stop_left = True
            elif len(deltaG_list[left_i-1:right_i+1][deltaG_list[left_i-1:right_i+1]<deltaG_thres]) > tolerance*(right_i-left_i+2):
                dormant_left = True
                if dormant_right==True or stop_right==True:
                    break
            else:
                dormant_left,dormant_right = False,False
                left_i -= 1
        if stop_right==False:
            if right_i+1>=len(st_list):
                stop_right = True
            elif st_list[right_i+1]-st_list[right_i]>max_gap:
                stop_right = True
            elif len(deltaG_list[left_i:right_i+2][deltaG_list[left_i:right_i+2]<deltaG_thres]) > tolerance*(right_i-left_i+2):
                dormant_right = True
                if dormant_left==True or stop_left==True:
                    break
            else:
                dormant_left,dormant_right = False,False
                right_i += 1

    
    if st_list[right_i] +1 - st_list[left_i] >= min_len:
    # if right_i-left_i >= min_len:
        # return pd.DataFrame({'ch':[ch],'st':[st_list[left_i]],'ed':[st_list[right_i] +1]})
        return [ch,st_list[left_i],st_list[right_i]+1]

    return None

def process_rows_parallel(ratio_df,deltaG_thres_seed,deltaG_thres, num_processes=None, Kd_column='Kd_log'):

    ratio_df_dict = ratio_split(ratio_df)
    footprint_df = pd.DataFrame({'ch':[],'st':[],'ed':[]})
    seed_index = list(ratio_df[ratio_df[Kd_column]>=deltaG_thres_seed].index)
    print('seed number', len(seed_index))

    partial_deNovoFootprinting = partial(deNovoFootprinting, ratio_df=ratio_df, ratio_df_dict=ratio_df_dict, deltaG_thres_seed=deltaG_thres_seed, deltaG_thres=deltaG_thres)

    # Determine the number of processes to use
    if num_processes is None:
        num_processes = multiprocessing.cpu_count()

    # Use multiprocessing.Pool with specified number of processes
    with multiprocessing.Pool(processes=num_processes) as pool:
        footprints = pool.map(partial_deNovoFootprinting, seed_index)

    footprints = [i for i in footprints if i!=None]
    footprints = sorted(footprints, key=lambda t: (t[0], t[1]))

    ch_list,st_list,ed_list = [],[],[]
    prev_ch_st_ed = ('chr1',0,0)
    for footprint in footprints:
        ch,st,ed = footprint
        ch_prev,st_prev,ed_prev = prev_ch_st_ed
        if ch==ch_prev and st<=ed_prev:
            continue
        prev_ch_st_ed = (ch,st,ed)
        ch_list.append(ch)
        st_list.append(st)
        ed_list.append(ed)

    footprint_df = pd.DataFrame({'ch':ch_list,'st':st_list,'ed':ed_list})
    footprint_df['st'] = np.array(footprint_df['st'],dtype=np.int32)
    footprint_df['ed'] = np.array(footprint_df['ed'],dtype=np.int32)

    return footprint_df


if __name__ == '__main__':
    in_f,out_f,deltaG_thres_seed,deltaG_thres,num_processes = sys.stdin.readline().strip().split(' ')
    print(in_f,out_f,deltaG_thres_seed,deltaG_thres,num_processes)
    deltaG_thres_seed,deltaG_thres = float(deltaG_thres_seed),float(deltaG_thres)
    num_processes = int(num_processes)
    ratio_df = pd.read_csv(in_f,sep='\t')
    ratio_df['Kd_log'] = -np.log10(ratio_df['Kd'])
    footprint_df = process_rows_parallel(ratio_df,deltaG_thres_seed,deltaG_thres, num_processes)
    footprint_df.to_csv(out_f,sep='\t',header=False,index=None)
    

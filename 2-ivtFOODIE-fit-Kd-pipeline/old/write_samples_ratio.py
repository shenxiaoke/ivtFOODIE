import tabix,pysam
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# import seaborn as sns
import pandas as pd
import glob,os

base_list = ['A','C','G','T','N']
compl_dict = {'A':'T','T':'A','C':'G','G':'C','N':'N'}

def complementary(context):
    return ''.join([compl_dict[i] for i in context])[::-1]

def filterDepth_siteAdjust(gz_dict, open_region, f_genome, min_depth=10):
    open_region_df = pd.read_csv(open_region, sep='\t',header=None)
    open_region_df.columns = ['ch','st','ed','_']
    
    ratio_dict = {}
    for index,region in open_region_df.iterrows():
        print(index, end='\r')
        ch,st,ed,_ = region
        # ref_seq = f_genome.fetch(ch,st,ed).upper()
        records_dict = {} 
        query_failed = False
        for gz_name in gz_dict:
            try:
                this_records = list(gz_dict[gz_name].query(ch,st,ed))
            except:
                query_failed = True
                break
            this_records_dict = {int(record[1]):record for record in this_records} # st: record
            records_dict[gz_name] = this_records_dict
        if query_failed == True:
            continue
        for pos in range(st,ed):
            this_ratio_dict = {}
            # base = ref_seq[pos-st]
            # if base == 'C' or base == 'G':
            depth_OK = True
            for gz_name in records_dict:
                if pos not in records_dict[gz_name]:
                    depth_OK = False
                    break
                con,uncon,ratio = int(records_dict[gz_name][pos][3]), int(records_dict[gz_name][pos][4]), float(records_dict[gz_name][pos][5])
                if con+uncon < min_depth:
                    depth_OK = False
                    break
                else:
                    this_ratio_dict[gz_name] = (ratio,con+uncon)
            if depth_OK == False:
                continue
            else:
                ratio_dict[(ch,pos)] = this_ratio_dict
        # print(ratio_dict)
        # break
    return ratio_dict


def write_ratio_f(ratio_dict,out_f,samples_from,samples_to):
    with open(out_f,'w') as fout:
        header = 'ch\tst\ted\t' + '\t'.join([f'ratio_{i}\tdepth_{i}' for i in range(samples_from,samples_to+1)]) + '\tratio_adjust\tdepth_adjust' + '\n'
        fout.write(header)
        for ch_pos in ratio_dict:
            if len(ch_pos)==2:
                ch,pos = ch_pos
                st,ed = pos,pos+1
            else:
                ch,st,ed = ch_pos
            line = f'{ch}\t{st}\t{ed}'
            for sample in range(samples_from,samples_to+1):
                ratio,depth = ratio_dict[ch_pos][sample]
                line += f'\t{ratio}\t{depth}'
            line += f'\t{ratio_dict[ch_pos]["adjust"][0]}\t{ratio_dict[ch_pos]["adjust"][1]}' + '\n'
            fout.write(line)


import argparse

genome_f='/dshare/xielab/analysisdata/TFC/shenk/important_data/chrom/genome/GRCh38/genome.fa'
f_genome = pysam.FastaFile(genome_f)

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description='create samples ratio table')

    # Add argument
    parser.add_argument('--TF', type=str, help='TF name')
    parser.add_argument('--samples_from', type=int, help='samples_from')
    parser.add_argument('--samples_to', type=int, help='samples_to')
    parser.add_argument('--sample_path', type=str, help='sample_path')
    parser.add_argument('--PC_path', type=str, help='PC_path')
    parser.add_argument('--out_path', type=str, help='out_path')
    parser.add_argument('--out_name', type=str, default=None, help='out_name')
    parser.add_argument('--open_region_f', type=str, help='open_region_f')
    parser.add_argument('--sample_index_position', type=int, default=-1, help='-1 means the last number of sample name represent the concentration, while -2 means the second to last')

    # Parse the arguments
    args = parser.parse_args()

    # Print the argument
    TF = args.TF
    samples_from = args.samples_from
    samples_to = args.samples_to
    sample_path = args.sample_path
    PC_path = args.PC_path
    out_path = args.out_path
    out_name = args.out_name
    open_region = args.open_region_f
    sample_index_position = args.sample_index_position

    gz_dict = {}

    gz_dict['adjust'] = tabix.open(PC_path)
    
    for i in range(samples_from,samples_to+1):
        for f in os.listdir(f'{sample_path}'):
            if str(i)!=f.split('-')[sample_index_position].split('_')[0]:
                continue
            print(f)
            for gz in os.listdir(f'{sample_path}/{f}'):
                if gz.split('.')[-1]!='gz':
                    continue
                gz_dict[i] = tabix.open(f'{sample_path}/{f}/{gz}')

    print(gz_dict)

    if len(gz_dict)!=1+(1+samples_to-samples_from):
        raise ValueError
    
    ratio_dict = filterDepth_siteAdjust(gz_dict, open_region, f_genome, min_depth=10)
    
    if out_name:
        out_f = f'{out_path}/{out_name}'
    else:
        out_f = f'{out_path}/{TF}.samples_ratio.csv'
    
    write_ratio_f(ratio_dict,out_f,samples_from,samples_to)

if __name__ == "__main__":
    main()
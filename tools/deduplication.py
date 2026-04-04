import numpy as np
import pysam,os,pickle,sys

bed_inf,bed_outf = sys.stdin.readline().strip().split(' ')

print(bed_inf)

def deduplication(bed_inf,bed_outf):
    # bed_inf is pair-in-one-line sorted (k1,1 k2,2n k3,3n) summary_score.bed
    bedfile = open(bed_inf,'r')
    bed_out = open(bed_outf,'w')
    prev_read_pair_pool = []
    prev_read_pair_chrom_st = ('chr1','0','0')
    in_pair_count = 0
    out_pair_count = 0
    # read pairs with the same positions
    pos_same = 0
    pos_same_read_different = 0
    for read in bedfile:
        in_pair_count += 1
        if in_pair_count % 100 == 0:
            print(in_pair_count,end='\r')
        # if in_pair_count % 1000000 == 0:
        #     break
        ch,st,ed,readname,pos,con,CorG,num,merge_st_ed = read.strip().split('\t')
        read_pair = [ch,st,ed,pos,con,CorG]
        read_pair_chrom_st = (ch,st,ed)
        if read_pair_chrom_st == prev_read_pair_chrom_st:
            duplication_mark = False
            for prev_read_pair in prev_read_pair_pool:
                if read_pair == prev_read_pair:
                    duplication_mark = True
                    break
            if duplication_mark==False:
                bed_out.write(read)
                prev_read_pair_pool.append(read_pair)
                out_pair_count += 1
                # how many read pairs are located in the same positions with different methylation
                pos_same_read_different += 1
            pos_same += 1
        else:
            bed_out.write(read)
            prev_read_pair_chrom_st = read_pair_chrom_st
            prev_read_pair_pool = [read_pair]
            out_pair_count += 1
    print(f'input read pairs:{in_pair_count},\t output read pairs:{out_pair_count}')
    print(f'pos_same_read_different:{pos_same_read_different},\t pos_same:{pos_same},\t ratio:{round(100*pos_same_read_different/pos_same,1)}%')
    bed_out.close()

deduplication(bed_inf,bed_outf)


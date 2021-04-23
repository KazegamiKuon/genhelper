import os
import zarr
import allel
import numpy as np
import time
import vcf_helper as vhelper
import ldhat_helper as ldhelper
from .config_genhelper import convert_values
import pandas as pd
from random import sample as sampling
from tqdm import tqdm
#Bio
from Bio import SeqIO as sqio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

def get_batch_snp_indexs(snp_indexs,size):
    data = snp_indexs.copy()
    batchs = []
    while len(data) > 0:
        batch = data[:size]
        batchs.append(batch)
        if size*0.9 > len(data):
            data = np.array([])
        else:
            data = data[int(size*0.9):]
    return batchs

def process_vcf(vcfpath,dataname,save_dir,bsize=96,bsnpsize=100000,scale=1e3):
    '''
    Input:
        vcfpath: đường dẫn file vcf
        dataname: tên dữ liệu/ tên file
        save_dir: thư mục lưu trữ file
        bsize: Số lượng sample mỗi batch
        bsnpsize: Số lượng snp mỗi sample mỗi batch
        scale: Chuyển đổi giá trị của locus xuống. Vd đang từ bp -> kb thì scale là 1e3
    '''
    zarrpath = vhelper.vcf_to_zarr(vcfpath)
    callset = zarr.open_group(zarrpath)
    #view information
    data_summary = vhelper.summary(vcfpath)
    snp_indexs, snp_poss = ldhelper.get_snp_index_snp_pos(callset)
    data_map, lhaplotypes, rhaplotypes = ldhelper.get_data_mapping(callset)

    # chia mẫu dựa trên sample
    nb_sample = data_summary[0]
    nb_snp = data_summary[4][0]

    sample_indexs = np.arange(nb_sample)
    batch_snp_indexs = get_batch_snp_indexs(np.arange(nb_snp),bsnpsize)
    batchs = ldhelper.get_sample_batch(sample_indexs,size=bsize)
    nb_batch = len(batchs)
    nb_batch_indexs = len(batch_snp_indexs)
    for i in tqdm(range(nb_batch)):
        for j in range(nb_batch_indexs):
            # Đường dẫn lưu file
            file_name = dataname+'_sbsize_{0}_nbsnp_{1}_b{2:03d}_p{3:03d}'.format(bsize,bsnpsize,i,j)
            sites_path = os.path.join(save_dir,file_name+'.sites')
            locs_path = os.path.join(save_dir,file_name+'.locs')
            ldhelper.create_sites_file(sites_path,data_map,lhaplotypes,rhaplotypes,batchs[i],snp_indexs[batch_snp_indexs[j]],callset.samples)
            ldhelper.create_locs_file(locs_path,snp_poss[batch_snp_indexs[j]],scale=scale)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcfpath', type=str,help='đường dẫn file vcf')
    parser.add_argument('-dataname', type=str,help='tên dữ liệu/ tên file dùng làm tên cho file output.')
    parser.add_argument('-save_dir', type=str,help='thư mục lưu trữ file output.')
    parser.add_argument('-bsize', type=int,default=96,help='Số lượng sample mỗi batch, default là 96')
    parser.add_argument('-bsnpsize', type=int,default=100000,help='Số lượng snp mỗi sample mỗi batch, default là 100000')
    parser.add_argument('-scale', type=int,default=convert_values._bp_to_kb,help='Chuyển đổi giá trị của locus xuống. Vd đang từ bp -> kb thì scale là 1e3. default là 1000)')
    args = parser.parse_args()
    process_vcf(args.vcfpath,args.dataname,args.save_dir,args.bsize,args.bsnpsize,args.scale)
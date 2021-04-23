import os
import zarr
import allel
import numpy as np
import time
import pandas as pd
from random import sample as sampling
from tqdm import tqdm
import sys

#Bio
from Bio import SeqIO as sqio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .config_genhelper import convert_values
# def get_range(poss: np.ndarray,length = 200):
#     tindexs =[]
#     to_bound_down = length*1000//20
#     to_bound_up = length*1000
#     bound_down = poss[0] - to_bound_down
#     bound_up =  bound_down +  to_bound_up
#     loop = (poss[-1] - bound_down)//to_bound_up # lấy só lần lặp
#     for i in range(1,loop+1):
#         cindexs = np.where((poss < bound_up) & (poss >= bound_down)) # Index của giá trị thỏa điều kiện
#         tindexs.append(dict({'pos_indexs':cindexs[0],'bound_down':bound_down,'bound_up':bound_up}))
#         bound_down = bound_up
#         bound_up = bound_down + to_bound_up
#     return loop, tindexs

def get_snp_index_snp_pos(database: zarr.hierarchy.Group):
    is_snp = database.variants.is_snp[:]
    indexs = np.argwhere(is_snp).flatten()
    poss = database.variants.POS[:][indexs]
    return indexs, poss

def get_data_mapping(database: zarr.hierarchy.Group): #,sample_indexs,snp_indexs
    stime = time.time()
    # lấy dữ liệu mapping
    data_map = np.vstack([database.variants.REF[:],database.variants.ALT[:].T]).T #[snp_indexs]
    lhaplotypes = database.calldata.GT[:,:,0]   # [snp_indexs][:,sample_indexs]
    rhaplotypes = database.calldata.GT[:,:,1]   # [snp_indexs][:,sample_indexs]
    print("--- %s seconds ---" % (time.time() - stime))
    return data_map, lhaplotypes, rhaplotypes

def create_sites_file(file_path,data_map,lhts,rhts,sample_indexs,snp_indexs,sample_names):
    nsamples = len(sample_indexs)
    hlength = len(snp_indexs)
    with open(file_path,mode='w') as seqs:
        seqs.write('{0} {1} {2}\n'.format(nsamples*2,hlength,1))
        for i in range(nsamples):
            sindex = sample_indexs[i] #index của sample i
            name = sample_names[sindex]
            
            # Haplotype trái của sample
            lhaplotype = data_map[snp_indexs,lhts.T[sindex][snp_indexs]]
            # Record của Haplotype trái
            lrecord = SeqRecord(seq=Seq(''.join(lhaplotype)),id='L'+name,name='L'+name,description='')
            # Lưu xuống file
            sqio.write(lrecord,handle=seqs,format='fasta')

            # Haplotype phải của sample
            rhaplotype = data_map[snp_indexs,rhts.T[sindex][snp_indexs]]
            # Record của Haplotype phải
            rrecord = SeqRecord(seq=Seq(''.join(rhaplotype)),id='R'+name,name='R'+name,description='')
            # Lưu xuống file
            sqio.write(rrecord,handle=seqs,format='fasta')

def create_sites_file_from_data(file_path,data):
    nsamples = len(data)
    hlength = len(data[0])
    with open(file_path,mode='w') as seqs:
        seqs.write('{0} {1} {2}\n'.format(nsamples,hlength,1))
        for i in range(nsamples):            
            # Record của Haplotype phải
            record = SeqRecord(seq=Seq(''.join(data[i])),id='sample_{0}'.format(i),name='sample_{0}'.format(i),description='')
            # Lưu xuống file
            sqio.write(record,handle=seqs,format='fasta')

def create_locs_file(file_path,poss,scale = convert_values._bp_to_kb):
    with open(file_path,mode='w') as locs:
        nloci = len(poss)
        locus = poss*scale
        length = np.round(locus[-1]+0.5,0)
        locs.write('{0} {1} L\n'.format(nloci,length))
        locs.write('\n'.join(locus.astype(str)))
    basename = os.path.basename(file_path)
    # truth_fpath = file_path.replace(basename,'truth_'+basename)
    # with open(truth_fpath,'w') as truthf:
    #     truthf.write(' '.join(poss.astype(str)))

def get_sample_batch(sample_indexs,size):
    # sample indexs, các giá trị không được trùng nhau
    samples = set(sample_indexs)
    sample_batch = []
    nb_loop = int(len(samples)/size)
    for i in range(nb_loop):
        sample = np.array(sampling(samples,size))
        sample_batch.append(sample)
        samples = samples.difference(sample)
    return sample_batch

if __name__ == '__main__':
    print('not supported')
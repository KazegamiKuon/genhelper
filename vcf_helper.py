import gzip
from operator import pos
import os
from types import FunctionType
import allel
import h5py
import zarr
# import numcodecs
import numpy as np
import pandas as pd
import time
# import fire
import argparse
from .config_class import vcf_zarr_config as vzconfig
import typing
from tqdm.notebook import tqdm

def get_df_exclude_par_region(df:pd.DataFrame,pos_col:str,pars:list)->pd.DataFrame:
    mask_array = None
    # get mask for par region
    for par in pars:
        temp = (df[pos_col]>= par[0])&(df[pos_col<=par[1]])
        if mask_array is None:
            mask_array= temp
        else:
            mask_array = mask_array & temp
        del temp
    # change mask to not par region
    mask_array = mask_array == False
    return df[mask_array].copy()

def consesus_genotype(data:tuple) -> str:
    # define special case, data only 2 elements
    temp_inter = set(data[0])
    list(map(lambda x: temp_inter.intersection_update(set(x)),data))
    nb_inter = len(temp_inter)
    if nb_inter >= 2:
        temp_inter = np.array(list(temp_inter))
    elif nb_inter == 1:
        temp_inter = np.array([-1,*list(temp_inter)])
    else:
        temp_inter = np.array([-1,-1])
    return '/'.join(list(map(vzconfig.gtmmap.get,temp_inter)))

def consesus_genotype_sample(data:list)->np.ndarray:
    temp = zip(*data)
    redata = np.array(list(map(consesus_genotype,temp)))
    return redata

def merge_df_variant_id(df_x:pd.DataFrame,df_y:pd.DataFrame):
    # id = 'id'
    # df_x[id] = df_x.apply(lambda row: ':'.join([row[vzconfig.chrom],str(row[vzconfig.position]),row[vzconfig.ref],row[vzconfig.alt]]),axis=1)
    # df_y[id] = df_y.apply(lambda row: ':'.join([row[vzconfig.chrom],str(row[vzconfig.position]),row[vzconfig.ref],row[vzconfig.alt]]),axis=1)
    return pd.merge(df_x,df_y,how='inner',on=[vzconfig.chrom,vzconfig.position,vzconfig.ref,vzconfig.alt])

def get_dataframe_variant_id(vtcallsets:typing.List[zarr.Group])-> pd.DataFrame:
    df_variant_id = None
    for i, vtcallset in enumerate(vtcallsets):
        tempdf = pd.DataFrame({
            vzconfig.chrom: vtcallset[vzconfig.chrom][:],
            vzconfig.position: vtcallset[vzconfig.position][:],
            vzconfig.ref: vtcallset[vzconfig.ref][:],
            vzconfig.alt: vtcallset[vzconfig.alt][:,0],            
            vzconfig.get_index_col(i): np.arange(vtcallset[vzconfig.chrom].shape[0])
        })
        if df_variant_id is None:
            df_variant_id = tempdf
        else:
            df_variant_id = merge_df_variant_id(df_variant_id,tempdf)
    return df_variant_id

def diploid_to_haploid_male(mapdata):
    male = mapdata[0]
    diploid = mapdata[1]
    temp = set(['0','1'])
    if male:
        temp.intersection_update(set(diploid.split('/')))
        if len(temp) == 0:
            return '.'
        else:
            return '/'.join(temp)
    return diploid

def consesus_genotype_vcf(vcf_files:typing.List[str],pars:list,vcf_consesus_prefix:str,check_male:FunctionType)->None:
    callsets = []
    data_idxs = []
    vcallsets = []
    for i,vcf_file in enumerate(tqdm(vcf_files,desc="create zarr callset")):
        zarr_path = vcf_to_zarr(vcf_file)
        callset = zarr.open_group(zarr_path)
        callsets.append(callset)
        vcallsets.append(callset.variants)
        data_idxs.append(vzconfig.get_index_col(i))
    inter_variantdf = get_dataframe_variant_id(vcallsets)
    ispar_mask_idxs = None
    for par in pars:
        positions = inter_variantdf[vzconfig.position].values
        mask_indexs = (( positions > par[0]) & (positions <= par[1]))
        if ispar_mask_idxs is not None:
            ispar_mask_idxs = ispar_mask_idxs | mask_indexs
        else:
            ispar_mask_idxs = mask_indexs
    inter_variantdf[vzconfig.par_col] = ispar_mask_idxs
    inter_variantdf.to_csv(vcf_consesus_prefix+'.csv.gz')
    gt_data = []
    for i, callset in enumerate(tqdm(callsets,desc="read gt data")):
        indexs = inter_variantdf[vzconfig.get_index_col(i)].values
        gt_data.append(callset.calldata[vzconfig.genotype].get_orthogonal_selection((indexs,slice(None),slice(None))))
    nb_variant = len(inter_variantdf)
    with gzip.open(vcf_consesus_prefix+'.vcf.gz','wt') as wf:
        sample_name = callsets[0].samples[:]
        males = list(map(check_male,sample_name))
        line = '\t'.join(sample_name)+'\n'
        wf.write(line)
        for i in tqdm(range(nb_variant),desc="intersection genotype"):
            samples = [gt[i] for gt in gt_data]
            items = consesus_genotype_sample(samples)
            if ispar_mask_idxs[i] == False:
                # zip data
                items = zip(males,items)
                # convert to haploid is sample is male
                items = list(map(diploid_to_haploid_male,items))
            line = '\t'.join(items)+'\n'
            wf.write(line)

def check_vcf_file(path):
    if not os.path.isfile(path):
        return False
    _, tfile = os.path.split(path)
    fname, fex1 = os.path.splitext(tfile)
    fname, fex2 = os.path.splitext(fname)
    if fex1 == '.gz' and (fex2 == '.vcf' or fex2 == '.bcf'):
        return True
    elif fex1 == '.vcf' or fex1 == '.bcf':
        return True
    return False

def count_alleles(variants : allel.VariantChunkedTable,ctype):
    assert ctype in ['snp','indel'] , "ctype must in ['snp','indel']"
    nb_0 = 0
    nb_1 = 0
    nb_2 = 0
    nb_var = variants.n_variants
    bsize = 10000
    loop = nb_var//bsize
    if nb_var%bsize != 0:
        loop += 1
    for i in range(loop):
        snps = variants.values.is_snp[bsize*i:bsize*i+bsize]
        nb_0 = nb_0 + np.sum(snps)
        if ctype == 'indel':
            temp = variants.values.altlen[bsize*i:bsize*i+bsize]
            nb_1 += np.sum(temp > 0)
            nb_2 += np.sum(temp < 0)
        if ctype == 'snp':
            refs = variants.values.REF[bsize*i:bsize*i+bsize][snps]
            alts = variants.values.ALT[bsize*i:bsize*i+bsize][snps]
            g = refs == 'G'
            a = refs == 'A'
            t = refs == 'T'
            c = refs == 'C'
            nb_1 += np.sum(alts[g,:] == 'A')
            nb_1 += np.sum(alts[a,:] == 'G')
            nb_1 += np.sum(alts[t,:] == 'C')
            nb_1 += np.sum(alts[c,:] == 'T')
            nb_2 += np.sum((alts[g,:] == 'T')+(alts[g,:] == 'C'))
            nb_2 += np.sum((alts[a,:] == 'T')+(alts[a,:] == 'C'))
            nb_2 += np.sum((alts[t,:] == 'A')+(alts[t,:] == 'G'))
            nb_2 += np.sum((alts[c,:] == 'A')+(alts[c,:] == 'G'))
    if ctype == 'indel':
        nb_0 = len(variants) - nb_0
    return nb_0, nb_1, nb_2

def get_zarr_path(path,in_zarr_folder=True):
    dpath, tfile = os.path.split(path)
    fname, fex = os.path.splitext(tfile)
    if in_zarr_folder:
        dpath, tfolder = os.path.split(dpath)
        dpath = os.path.join(dpath,'zarr')
    return os.path.join(dpath,fname+'.zarr')

def vcf_to_zarr(vcfpath,in_zarr_folder=True):
    zarrpath = get_zarr_path(vcfpath,in_zarr_folder)
    if not os.path.isdir(zarrpath):
        start_time = time.time()
        allel.vcf_to_zarr(vcfpath, zarrpath , fields='*', overwrite=True)
        print("--- %s seconds ---" % (time.time() - start_time))
    return zarrpath

def summary(vcfpath):
    # vcfpath = '1KGP.chr1.10M.hg38_phased.vcf.gz'
    zarrpath = ''
    if os.path.isfile(vcfpath):
        zarrpath = vcf_to_zarr(vcfpath)
    else:
        zarrpath = get_zarr_path(vcfpath)
    if not os.path.isdir(zarrpath):
        print("Zarr data of vcf is undifined. Phease check vcf path is True or Zarr already haved.")
        return
    
    start_time = time.time()

    callset = zarr.open_group(zarrpath)
    variants = allel.VariantChunkedTable(callset.variants)
    nb_sample = callset.samples.size
    nb_chormosome = len(np.unique(variants.values.CHROM))
    nb_variant = nb_record = variants.n_variants
    print('No. samples: {0}'.format(nb_sample))
    print('No. chromosome: ',nb_chormosome)
    print('No. variant: {0}'.format(nb_variant))
    print('No. record: {0}'.format(nb_record))

    nb_snp, nb_trnt, nb_trnv = count_alleles(variants,'snp')
    nb_indel, nb_ins, nb_del = count_alleles(variants,'indel')
    print('\n')
    print('No. SNP: {0}  [{1}/{2}] - [transitions/tranvertions]'.format(nb_snp,nb_trnt,nb_trnv))
    print('No. INDEL: {0}  [{1}/{2}] - [insertions/deletions]'.format(nb_indel,nb_ins,nb_del))
    print('\n')
    print("--- %s seconds ---" % (time.time() - start_time))
    return nb_sample, nb_chormosome, nb_variant, nb_record, [nb_snp,nb_trnt,nb_trnv], [nb_indel,nb_ins,nb_del]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # parser.add_argument('--size', type=int, default=256)
    # parser.add_argument('-type', type=str)
    # parser.add_argument('--name', type=str)
    parser.add_argument('-vcfpath', type=str)

    args = parser.parse_args()
    data = summary(args.vcfpath)
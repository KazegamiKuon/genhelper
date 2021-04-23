import os
import allel
import h5py
import zarr
import numcodecs
import numpy as np
import time
import fire
import argparse

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

def get_zarr_path(path):
    dpath, tfile = os.path.split(path)
    fname, fex = os.path.splitext(tfile)
    dpath, tfolder = os.path.split(dpath)
    dpath = os.path.join(dpath,'zarr')
    return os.path.join(dpath,fname+'.zarr')

def vcf_to_zarr(vcfpath):
    zarrpath = get_zarr_path(vcfpath)
    if not os.path.isdir(zarrpath):
        start_time = time.time()
        allel.vcf_to_zarr(vcfpath, zarrpath , fields='*', overwrite=True)
        print("--- %s seconds ---" % (time.time() - start_time))
    return zarrpath

def summary(vcfpath):
    # vcfpath = '1KGP.chr1.10M.hg38_phased.vcf.gz'
    zarrpath = vcf_to_zarr(vcfpath)
    
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
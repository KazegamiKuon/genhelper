class VcfToZarrConfig:
    def __init__(self) -> None:
        self.chrom = "CHROM"
        self.position = "POS"
        self.ref = "REF"
        self.alt = "ALT"
        self.par_col = "ispar"
        self.genotype = "GT"
        self.gtmmap = dict()
        self.gtmmap[-1] = '.'
        self.gtmmap['.'] = -1
        self.gtmmap[0] = '0'
        self.gtmmap['0'] = 0
        self.gtmmap[1] = '1'
        self.gtmmap['1'] = 1
        pass

    def get_index_col(self,i:int)->str:
        return 'index_{}'.format(i)

vcf_zarr_config = VcfToZarrConfig()
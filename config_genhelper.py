import numpy as np

class ConvertValues():
    def __init__(self):
        self._4Ne = 4*10000
        self._M_to_cM = 1/0.01
        self._4Ner_to_cM = (1/self._4Ne)*self._M_to_cM
        self._bp_to_kb = 1e-3
        self._bp_to_Mb = 1e-6
        self._kb_to_Mb = 1e-3
        self._Mb_to_kb = 1e3
        self._4Ner_per_kb_to_cM_per_Mb = (1/self._4Ne)*self._M_to_cM/self._kb_to_Mb
        self.convert_names = []

convert_values = ConvertValues()
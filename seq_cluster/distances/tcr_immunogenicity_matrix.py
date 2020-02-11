import math
import ntpath
import os
from string import ascii_uppercase

import pandas as pd
from numpy import mean

from seq_cluster.distances.distance_matrix import DistanceMatrix


class TCRImmunogenicityMatrix(DistanceMatrix):
    def __init__(self, csv_data):
        super(TCRImmunogenicityMatrix, self).__init__()
        self._init_data(csv_data)

    def _init_data(self, csv_data):
        for c0 in ascii_uppercase:
            for c1 in ascii_uppercase:
                key = c0 + c1
                if key in csv_data.index:
                    val = csv_data[key]
                    self._data[key] = val

    def get_distance(self, aa1, aa2):
        out = self._data.get((aa1 + aa2), None)
        if out is None:
            out = self._data.get((aa2 + aa1), None)
        return out

    def get_stats(self):
        out = ''
        identical_val=[]
        for id, value in self._data.iteritems():
            c0 = id[0:1]
            c1 = id[1:2]
            if c0 == c1:
                identical_val.append(value)
        pg = self.get_distance('P','G')
        max_v = max(self._data.values())
        max_k = self._get_key_for_value(max_v)
        min_v = min(self._data.values())
        min_k = self._get_key_for_value(min_v)
        out = 'Avg. identicals: ' + str(mean(identical_val)) + \
              ', max:' + max_k + '_' + str(max_v)  + \
              ', min:' + min_k + '_' + str(min_v) + ', PG: ' + str(pg)
        return out

    def _get_key_for_value(self, value):
        for id, val in self._data.iteritems():
            if value == val: return id
        return None

class TCRImmunogenicityMatrices(object):
    # From: Quantitative Prediction of the Landscape of T Cell Epitope Immunogenicity in Sequence Space, Masato Ogishi and Hiroshi Yotsuyanagi
    _MATRICES_DIRNAME = 'matrices'

    def __init__(self):
        self._items = {}
        fn = TCRImmunogenicityMatrices.get_matrix_data()
        df = pd.read_csv(fn, sep=',')
        for row in df.iterrows():
            tmp, data = row
            m = TCRImmunogenicityMatrix(data)
            self._items[data['AAIndexID']] = m

    @staticmethod
    def get_matrix_data():
        out = os.path.realpath(__file__)
        out = ntpath.dirname(out)
        out = os.path.join(out, TCRImmunogenicityMatrices._MATRICES_DIRNAME, 'tcr_immunogenicity.csv')
        return out

    def get_keys(self):
        return self._items.keys()

    def get_matrix(self, name):
        return self._items.get(name)


ms  =TCRImmunogenicityMatrices()
for n in ms.get_keys():
    m = ms.get_matrix(n)
    s = m.get_stats()
    print s
print 2
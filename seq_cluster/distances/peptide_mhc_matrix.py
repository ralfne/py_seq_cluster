import ntpath
import os

import pandas as pd
import numpy as np

from seq_cluster.distances.distance_matrix import DistanceMatrix


class PeptideMhcMatrix(DistanceMatrix):
    _MATRICES_DIRNAME = 'matrices'

    def __init__(self):
        super(PeptideMhcMatrix, self).__init__()
        df = self._load_matrix()
        df = self._transform_to_distances(df)
        self._data = self._convert_to_dict(df)

    def _load_matrix(self):
        fn = PeptideMhcMatrix.get_matrix_data()
        df = pd.read_csv(fn, sep='\t', index_col=0)
        return df

    def _transform_to_distances(self, df):
        max_value = df.max().max()
        df = df * -1
        df = df + max_value
        return df

    def _convert_to_dict(self, df):
        out = {}
        df = df.where(np.triu(np.ones(df.shape)).astype(np.bool))
        df = df.stack().reset_index()
        df.columns = ['Row', 'Column', 'Value']
        for key, row in df.iterrows():
            c0 = row['Row']
            c1 = row['Column']
            v = row['Value']
            out[c0, c1] = v
        return out

    def get_distance(self, aa1, aa2):
        out = self._data.get((aa1, aa2), None)
        if out is None:
            out = self._data.get((aa2, aa1), None)
        return out

    @staticmethod
    def get_matrix_data():
        out = os.path.realpath(__file__)
        out = ntpath.dirname(out)
        out = os.path.join(out, PeptideMhcMatrix._MATRICES_DIRNAME, 'peptide_mhc_matrix.txt')
        return out

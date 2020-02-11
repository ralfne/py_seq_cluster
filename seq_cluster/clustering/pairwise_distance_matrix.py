from abc import ABCMeta, abstractmethod
from copy import deepcopy
import Bio.SubsMat.MatrixInfo as biomat
from similarity.levenshtein import Levenshtein
import pandas as pd
import numpy as np
from similarity.weighted_levenshtein import WeightedLevenshtein
from similarity.weighted_levenshtein import CharacterSubstitutionInterface


class PairwiseDistanceMatrix(object):
    def __init__(self, pairwise_comparisons, calculator):
        self._calculator = calculator
        self._pairwise_comparisons = pairwise_comparisons
        self._calculate_distances()

    def _calculate_distances(self):
        for pair in self._pairwise_comparisons:
            self._calculate_distances_for_pair(pair)

    def _calculate_distances_for_pair(self, pair):
        item1 = pair.get_item1()
        item2 = pair.get_item2()
        dist = self._calculator.calculate(item1, item2)
        pair.set_distance(dist)

    def get_as_dataframe(self):
        names = self._get_unique_names()
        data = []
        for row_name in names:
            row = []
            for col_name in names:
                if row_name  == col_name:
                    v = 0
                else:
                    v = self._get_distance(row_name, col_name)
                row.append(v)
            data.append(row)
        out = pd.DataFrame(data=data, index=deepcopy(names), columns=deepcopy(names))
        return out

    def _get_distance(self, name1, name2):
        for p in self._pairwise_comparisons:
            if p.matches(name1, name2, respect_ordering=False):
                return p.get_distance()
        return None

    def _get_unique_names(self):
        item1 = self._pairwise_comparisons[0].get_item1()
        out = [item1.name]
        for p in self._pairwise_comparisons:
            item1 = p.get_item1()
            item2 = p.get_item2()
            if item1.name == out[0]:
                out.append(item2.name)
        return out

    @staticmethod
    def combine(lower_triangle_distance_dataframe, upper_triangle_distance_dataframe):
        X_upper = upper_triangle_distance_dataframe.values
        v_upper = X_upper[np.triu_indices(X_upper.shape[0], k=0)]
        X_lower = lower_triangle_distance_dataframe.values
        v_lower = X_lower[np.tril_indices(X_lower.shape[0], k=0)]
        size_X = len(lower_triangle_distance_dataframe)
        out = np.zeros((size_X, size_X))
        out[np.triu_indices(out.shape[0], k=0)] = v_upper
        out[np.tril_indices(out.shape[0], k=0)] = v_lower
        out = pd.DataFrame(out, index=lower_triangle_distance_dataframe.index, columns=upper_triangle_distance_dataframe.columns)
        return out

from abc import ABCMeta, abstractmethod


class DistanceMatrix(object):
    __metaclass__ = ABCMeta

    def __init__(self):
        self._data = {}

    @abstractmethod
    def get_distance(self, aa1, aa2): pass

    def get_min_max_matrix_values(self):
        min_value = min(self._data.values())
        max_value = max(self._data.values())
        return min_value, max_value
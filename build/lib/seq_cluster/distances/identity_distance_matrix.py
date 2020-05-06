from string import ascii_uppercase
from seq_cluster.distances.distance_matrix import DistanceMatrix


class IdentityDistanceMatrix(DistanceMatrix):
    def __init__(self):
        super(IdentityDistanceMatrix, self).__init__()
        self._data = self._get_all_combinations()

    def _get_all_combinations(self):
        out = {}
        for c0 in ascii_uppercase:
            for c1 in ascii_uppercase:
                key = c0, c1
                if c0 == c1:
                    out[key] = 0
                else:
                    out[key] = 1
        return out

    def get_distance(self, aa1, aa2):
        key = (aa1, aa2)
        out = self._data.get(key, None)
        if out is None: raise ValueError('Key not found: ' + str(key))
        return out
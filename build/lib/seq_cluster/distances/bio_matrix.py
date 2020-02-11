from seq_cluster.distances.distance_matrix import DistanceMatrix
import Bio.SubsMat.MatrixInfo as biomat


class BioMatrix(DistanceMatrix):
    def __init__(self, name):
        super(BioMatrix, self).__init__()
        data = getattr(biomat, name)
        self._data = self._transform_matrix(data)
        #substitution_matrix = biomat.blosum62.values()

    def _transform_matrix(self, data):
        out = {}
        max_value = max(data.values())
        min_value = min(data.values())
        for key, value in data.iteritems():
            v = value * -1
            v = v + max_value
            out[key] = v
        return out

    def get_distance(self, aa1, aa2):
        out = self._data.get((aa1, aa2), None)
        if out is None:
            out = self._data.get((aa2, aa1), None)
        return out

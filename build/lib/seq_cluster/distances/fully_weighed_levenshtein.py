# Copyright (c) 2018 luozhouyang
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
from copy import deepcopy

from logger.Logger import StdOutLogger

from seq_cluster.distances.peptide_mhc_matrix import PeptideMhcMatrix


class CharacterInsDelInterface:

    def deletion_cost(self, c):
        raise NotImplementedError()

    def insertion_cost(self, c):
        raise NotImplementedError()

class CharacterInsDel:
    def __init__(self, deletion_cost, insertion_cost):
        self._deletion_cost = deletion_cost
        self._insertion_cost = insertion_cost

    def deletion_cost(self):
        return self._deletion_cost

    def insertion_cost(self):
        return self._insertion_cost


class FullyWeightedLevenshtein(object):
    def __init__(self, distance_matrix, character_ins_del=None, logger=StdOutLogger(verbose=False)):
        self.character_ins_del = character_ins_del
        self.distance_matrix  = distance_matrix
        self._logger = logger
        self._cache_v0 = False
        self._v0s = {}

    def distance(self, s0, s1):
        if s0 is None:
            raise TypeError("Argument s0 is NoneType.")
        if s1 is None:
            raise TypeError("Argument s1 is NoneType.")
        if s0 == s1:
            return 0.0
        if len(s0) == 0:
            return len(s1)
        if len(s1) == 0:
            return len(s0)

        if self._cache_v0:
            v0, v1 = self._get_v0_v1(s1)
        else:
            v0, v1 = [0.0] * (len(s1) + 1), [0.0] * (len(s1) + 1)
            v0[0] = 0
            for i in range(1, len(v0)):
                v0[i] = v0[i - 1] + self._insertion_cost()
                # v0[i] = v0[i - 1] + self._insertion_cost(s1[i - 1])

        for i in range(len(s0)):
            s1i = s0[i]
            deletion_cost = self._deletion_cost()
            # deletion_cost = self._deletion_cost(s1i)
            v1[0] = v0[0] + deletion_cost

            for j in range(len(s1)):
                s2j = s1[j]
                cost = 0
                if s1i != s2j:
                    cost = self.distance_matrix.get_distance(s1i, s2j)
                # insertion_cost = self._insertion_cost(s2j)
                insertion_cost = self._insertion_cost()
                v1[j + 1] = min(v1[j] + insertion_cost, v0[j + 1] + deletion_cost, v0[j] + cost)
            v0, v1 = v1, v0

        return v0[len(s1)]

    def _get_v0_v1(self, s1):
        v0_v1 = self._v0s.get(s1, None)
        if v0_v1 is None:
            v0_v1 = self._create_v0_v1(s1)
            self._v0s[s1] = v0_v1
        return deepcopy(v0_v1)

    def _create_v0_v1(self, s1):
        v0, v1 = [0.0] * (len(s1) + 1), [0.0] * (len(s1) + 1)
        v0[0] = 0
        for i in range(1, len(v0)):
            v0[i] = v0[i - 1] + self._insertion_cost()
        return v0, v1

    def _insertion_cost(self):
        if self.character_ins_del is None:
            return 1.0
        return self.character_ins_del.insertion_cost()

    def _deletion_cost(self):
        if self.character_ins_del is None:
            return 1.0
        return self.character_ins_del.deletion_cost()



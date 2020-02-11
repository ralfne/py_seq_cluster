from abc import ABCMeta, abstractmethod
from copy import deepcopy

from Logger import StdOutLogger
from similarity.levenshtein import Levenshtein
import pandas as pd
import numpy as np
from similarity.weighted_levenshtein import WeightedLevenshtein
from similarity.weighted_levenshtein import CharacterSubstitutionInterface

from seq_cluster.distances.fully_weighed_levenshtein import FullyWeightedLevenshtein


# class CharacterSubstitution(CharacterSubstitutionInterface):
#     def __init__(self, substitution_matrix):
#         self._substitution_matrix = substitution_matrix
#         self._max=max(substitution_matrix.values())
#
#     def cost(self, c0, c1):
#         key = c0, c1
#         out = self._substitution_matrix.get(key, None)
#         if out is None:
#             key = c1, c0
#             out = self._substitution_matrix.get(key, None)
#         #print (out)
#         out = out - self._max
#         #print (out*-1)
#         return out*-1
#
#
# class Insdel(object):
#     def __init__(self, insertion=1, deletion=1):
#         self._insertion = insertion
#         self._deletion = deletion
#
#     def insertion_cost(self, c):
#         return self._insertion
#
#     def deletion_cost(self, c):
#         return self._deletion


class SequenceDistanceCalculator():
    def __init__(self,  levenshtein, v_j_cost_calculator, logger=StdOutLogger(verbose=False)):
        self._levenshtein = levenshtein
        self._v_j_cost_calculator = v_j_cost_calculator
        self._logger = logger

    def calculate(self, comparison_item1, comparison_item2):
        msg = 'Calculating distance for %s vs. %s ' %(comparison_item1.name, comparison_item2.name)
        seq1 = comparison_item1.sequence
        seq2 = comparison_item2.sequence
        out = self._levenshtein.distance(seq1, seq2)
        msg += '; Levensthein: %s' % (str(out))
        vj_cost = self._v_j_cost_calculator.calculate(comparison_item1, comparison_item2)
        msg += '; v_j_cost: %s' % (str(vj_cost))
        out += vj_cost
        self._logger.log(msg, onlyIfVerbose=True)
        return out


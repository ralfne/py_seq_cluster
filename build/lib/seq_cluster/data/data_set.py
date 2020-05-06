import uuid
from copy import deepcopy
import immune_receptor_utils.genes as irg
import pandas as pd


class ComparisonItem(object):
    def __init__(self, name, seq1=None, seq2=None, v_gene=None, j_gene=None, aux=None):
        self.name = name
        self.sequence1 = seq1
        self.sequence2 = seq2
        self.v_gene = v_gene
        self.j_gene = j_gene
        self.aux = aux

    def set_as_aux(self, aux):
        self.aux = aux
        self.sequence = None
        self.v_gene = None
        self.j_gene = None


class PairwiseComparison(object):
    def __init__(self, item1, item2):
        self._item1 = item1
        self._item2 = item2
        self._distance = None

    def set_distance(self, distance):
        self._distance = distance

    def get_distance(self):
        return self._distance

    def get_item1(self):
        return self._item1

    def get_item2(self):
        return self._item2

    def matches(self, name1, name2, respect_ordering=False):
        if (self._item1.name == name1) and (self._item2.name == name2):
            return True
        if not respect_ordering:
            if (self._item2.name == name1) and (self._item1.name == name2 ):
                return True
        return False


class Dataset(object):
    def __init__(self, filename, sequence1_col, sequence2_col, name_cols=None, v_genes_col=None, j_genes_col= None, aux_col=None):
        self._data = pd.read_csv(filename, sep='\t')
        self._sequence1_col = sequence1_col
        self._sequence1_col_index = None
        self._sequence2_col = sequence2_col
        if self._sequence2_col == '': self._sequence2_col = None
        self._sequence2_col_index = None
        self._name_cols = name_cols
        self._aux_col = aux_col
        self._v_genes_col = v_genes_col
        self._j_genes_col = j_genes_col
        self._assert_name_cols_are_unique()

    def _assert_name_cols_are_unique(self):
        if self._name_cols is not None:
            tmp = {}
            for i in range(len(self._data)):
                name = self.get_name(i)
                tmp[name] = i
            if len(tmp) != len(self._data):
                raise ValueError('Error: Name columns do not define unique row names.')

    def get_sequences1(self, row_index):
        out = self._data.iloc[row_index, self._get_sequence1_col_index()]
        return out

    def get_sequences2(self, row_index):
        if self._sequence2_col is None: return None
        out = self._data.iloc[row_index, self._get_sequence2_col_index()]
        return out

    def _get_sequence1_col_index(self):
        if self._sequence1_col_index is  None:
            self._sequence1_col_index = self._data.columns.get_loc(self._sequence1_col)
        return self._sequence1_col_index

    def _get_sequence2_col_index(self):
        if self._sequence2_col_index is  None:
            self._sequence2_col_index = self._data.columns.get_loc(self._sequence2_col)
        return self._sequence2_col_index

    def get_aux(self, row_index):
        row = self._data.iloc[row_index]
        out = row[self._aux_col]
        return out

    def get_v_gene(self, row_index):
        out = None
        if self._v_genes_col is not None:
            row = self._data.iloc[row_index]
            s = row[self._v_genes_col]
            out = irg.GenePossibilities(s, False)
        return out

    def get_j_gene(self, row_index):
        out = None
        if self._j_genes_col is not None:
            row = self._data.iloc[row_index]
            s = row[self._j_genes_col]
            out = irg.GenePossibilities(s, False)
        return out

    def get_name(self, row_index):
        if self._name_cols is None:
            out = str(row_index)
        else:
            row = self._data.iloc[row_index]
            out = ''
            for n in self._name_cols:
                if len(out)>0: out += '_'
                out += str(row[n])
        return out

    def get_pairwise(self):
        seq_out = []
        aux_out = []
        if self._aux_col is None:
            aux_out = None
        for i in range(len(self._data)-1):
            for j in range(i+1, len(self._data)):
                item1 = self._create_sequence_comparison_item(i)
                item2 = self._create_sequence_comparison_item(j)
                pwc = PairwiseComparison(item1, item2)
                seq_out.append(pwc)
                if self._aux_col is not None:
                    aux1 = self.get_aux(i)
                    aux2 = self.get_aux(j)
                    item1_aux = deepcopy(item1)
                    item1_aux.set_as_aux(aux1)
                    item2_aux = deepcopy(item2)
                    item2_aux.set_as_aux(aux2)
                    pwc_aux = PairwiseComparison(item1_aux, item2_aux)
                    aux_out.append(pwc_aux)
        return seq_out, aux_out

    def _create_sequence_comparison_item(self, row_index):
        seq1 = self.get_sequences1(row_index)
        seq2 = self.get_sequences2(row_index)
        name = self.get_name(row_index)
        v_gene = self.get_v_gene(row_index)
        j_gene = self.get_j_gene(row_index)
        out = ComparisonItem(name, seq1, seq2, v_gene, j_gene, None)
        return out

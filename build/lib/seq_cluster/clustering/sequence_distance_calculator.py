from logger.Logger import StdOutLogger


class SequenceDistanceCalculator():
    def __init__(self,  levenshtein, v_j_cost_calculator, logger=StdOutLogger(verbose=False)):
        self._levenshtein = levenshtein
        self._v_j_cost_calculator = v_j_cost_calculator
        self._logger = logger

    def calculate(self, comparison_item1, comparison_item2):
        msg = 'Calculating distance for %s vs. %s ' %(comparison_item1.name, comparison_item2.name)
        item1_seq1 = comparison_item1.sequence1
        item1_seq2 = comparison_item1.sequence2
        item2_seq1 = comparison_item2.sequence1
        item2_seq2 = comparison_item2.sequence2
        out = self._levenshtein.distance(item1_seq1, item2_seq1)
        if item1_seq2 is not None:
            out2 = self._levenshtein.distance(item1_seq2, item2_seq2)
            out += out2
        msg += '; Levensthein: %s' % (str(out))
        vj_cost = self._v_j_cost_calculator.calculate(comparison_item1, comparison_item2)
        msg += '; v_j_cost: %s' % (str(vj_cost))
        out += vj_cost
        self._logger.log(msg, onlyIfVerbose=True)
        return out


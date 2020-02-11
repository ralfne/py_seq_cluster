class SimpleAbsDistanceCalculator(object):
    def calculate(self, comparison_item1, comparison_item2):
        out = abs(comparison_item1.aux - comparison_item2.aux)
        return out
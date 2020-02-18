from Levenshtein.StringMatcher import StringMatcher


class FastLevenshteinCalculator(object):

    def distance(self, seq1, seq2):
        matcher = StringMatcher(seq1=seq1, seq2=seq2)
        out = matcher.distance()
        return out

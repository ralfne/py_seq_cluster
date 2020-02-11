class VJCostCalculator(object):
    def __init__(self, unequal_v_gene_cost=15, unequal_j_gene_cost=10):
        self._unequal_v_gene_cost = unequal_v_gene_cost
        self._unequal_j_gene_cost = unequal_j_gene_cost

    def calculate(self, comparison_item1, comparison_item2):
        out = 0
        if comparison_item1.v_gene is not None:
            if comparison_item2.v_gene is None: raise ValueError('Unequal V gene status')
            if not comparison_item1.v_gene.overlap_exists(comparison_item2.v_gene, include_allele=False):
                out += self._unequal_v_gene_cost
        if comparison_item1.j_gene is not None:
            if comparison_item2.j_gene is None: raise ValueError('Unequal J gene status')
            if not comparison_item1.j_gene.overlap_exists(comparison_item2.j_gene, include_allele=False):
                out += self._unequal_j_gene_cost
        return out
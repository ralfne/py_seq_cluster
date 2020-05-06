from seq_cluster.distances.bio_matrix import BioMatrix
from seq_cluster.distances.fully_weighed_levenshtein import FullyWeightedLevenshtein, CharacterInsDel
from seq_cluster.distances.identity_distance_matrix import IdentityDistanceMatrix
from seq_cluster.clustering.pairwise_distance_matrix import PairwiseDistanceMatrix
from seq_cluster.clustering.simple_abs_distance_calculator import SimpleAbsDistanceCalculator
from seq_cluster.clustering.v_j_cost_calculator import VJCostCalculator
from seq_cluster.clustering.sequence_distance_calculator import SequenceDistanceCalculator
from seq_cluster.data.data_set import Dataset
from seq_cluster.distances.fast_levenshtein_calculator import FastLevenshteinCalculator
from seq_cluster.distances.peptide_mhc_matrix import PeptideMhcMatrix
from seq_cluster.visualization.clustermap_visualizer import ClustermapVisualizer


def run(args, logger):
    logger.log('Loading sequences...', includeTimestamp=True, onlyIfVerbose=False)
    ds = Dataset(args.infile, sequence1_col=args.sequence1_col, sequence2_col=args.sequence2_col, name_cols=args.name,
            v_genes_col=args.v_gene, j_genes_col=args.j_gene, aux_col=args.aux)
    logger.log('Initializing pairs...', includeTimestamp=True, onlyIfVerbose=False)
    seq_pairs, aux_pairs = ds.get_pairwise()
    v_j_cost_calculator = VJCostCalculator(unequal_v_gene_cost=args.v_gene_cost,
                                            unequal_j_gene_cost=args.j_gene_cost)

    if args.subst_matrix == 'identity_matrix':
        subst_matrix = IdentityDistanceMatrix()
    elif args.subst_matrix == 'peptide_mhc_matrix':
        subst_matrix = PeptideMhcMatrix()
    elif args.subst_matrix == 'fast_levenshtein':
        subst_matrix = None
    else:
        subst_matrix = BioMatrix(args.subst_matrix)

    # min_v, max_v = subst_matrix.get_min_max_matrix_values()
    # msg = 'Using matrix: %s; min and max values: %s and %s' % (args.subst_matrix, min_v, max_v)
    # logger.log(msg,onlyIfVerbose=False)
    #
    # insdel = CharacterInsDel(args.del_cost, args.ins_cost)
    # levenshtein = FullyWeightedLevenshtein(subst_matrix, insdel, logger)
    # seq_calc = SequenceDistanceCalculator(levenshtein, v_j_cost_calculator, logger)
    if subst_matrix is None:
        seq_calc = _create_fast_levenshtein_calculator(v_j_cost_calculator, logger)
    else:
        seq_calc = _create_default_levenshtein_calculator(args, subst_matrix, v_j_cost_calculator, logger)

    logger.log('Calculating pairwise distances...', includeTimestamp=True, onlyIfVerbose=False)
    seq_matrix = PairwiseDistanceMatrix(seq_pairs, seq_calc)
    logger.log('Formatting pairwise distances...', includeTimestamp=True, onlyIfVerbose=False)
    seq_df = seq_matrix.get_as_dataframe()

    aux_df = None
    if aux_pairs is not None:
        aux_calc = SimpleAbsDistanceCalculator()
        aux_matrix = PairwiseDistanceMatrix(aux_pairs, aux_calc)
        aux_df = aux_matrix.get_as_dataframe()

    logger.log('Running clustermap...', includeTimestamp=True, onlyIfVerbose=False)

    fig_sizes = args.fig_size.split(',')

    viz = ClustermapVisualizer(seq_df, args.seq_linkage, aux_df, args.aux_linkage, scaling=0.2,
                               fig_size_x=int(fig_sizes[0].strip()), fig_size_y=int(fig_sizes[1].strip()),
                                fig_size_font=int(fig_sizes[2].strip()))
    viz.run(title=None, filename=args.outfile)
    if args.outfile is None:
        import matplotlib.pyplot as plt
        plt.show()

def _create_default_levenshtein_calculator(args, subst_matrix, v_j_cost_calculator, logger):
    min_v, max_v = subst_matrix.get_min_max_matrix_values()
    msg = 'Using matrix: %s; min and max values: %s and %s' % (args.subst_matrix, min_v, max_v)
    logger.log(msg,onlyIfVerbose=False)

    insdel = CharacterInsDel(args.del_cost, args.ins_cost)
    levenshtein = FullyWeightedLevenshtein(subst_matrix, insdel, logger)
    seq_calc = SequenceDistanceCalculator(levenshtein, v_j_cost_calculator, logger)
    return seq_calc

def _create_fast_levenshtein_calculator(v_j_cost_calculator, logger):
    levenshtein = FastLevenshteinCalculator()
    seq_calc = SequenceDistanceCalculator(levenshtein, v_j_cost_calculator, logger)
    return seq_calc

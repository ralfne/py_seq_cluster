import argparse
import Bio.SubsMat.MatrixInfo as biomat

def create_arguments():
    parser = argparse.ArgumentParser(description='Sequence clustering',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('infile', type=str, metavar='',
                        help='Input filename of the tab-separated data file; '
                             'use "show matrices" to display all available matrices.')
    parser.add_argument('sequence1_col', type=str, metavar='',
                        help='Column name for the sequence_1 clustering parameter')
    parser.add_argument('sequence2_col', type=str, metavar='', default='',
                        help='Column name for the sequence_2 clustering parameter (enter "" if not used)')
    parser.add_argument('-n', '--name', action='append', metavar='',
                        help='Column name specifying the item name. Multiple -n options are supported; '
                             'the resulting values will be concatenated')
    parser.add_argument('-x', '--aux', default=None, type=str, metavar='',
                        help='Column name for the auxillary clustering parameter')
    parser.add_argument('-s', '--seq_linkage', default='average', type=str, metavar='',
                        help='Linkage method for sequence clustering')
    parser.add_argument('-a', '--aux_linkage', default='average', type=str, metavar='',
                        help='Linkage method for auxillary clustering')
    parser.add_argument('-v', '--v_gene', default=None, type=str, metavar='',
                        help='Column name for V genes')
    parser.add_argument('-y', '--v_gene_cost', default=None, type=float, metavar='',
                        help='Cost associated with unequal V genes')
    parser.add_argument('-j', '--j_gene', default=None, type=str, metavar='',
                        help='Column name for J genes')
    parser.add_argument('-i', '--j_gene_cost', default=None, type=float, metavar='',
                        help='Cost associated with unequal J genes')
    parser.add_argument('-o', '--outfile', default=None, type=str, metavar='',
                        help='Output filename of the tab-separated data file')
    parser.add_argument('-m', '--subst_matrix', default='peptide_mhc_matrix', type=str,  metavar='',
                        help='Substitution matrix for calculating sequence distances, '
                             'Where applicable, the matrix will be transformed from a substitution matrix '
                            '(higher values=more similar) into a distance matrix (higher values=less similar)')
    parser.add_argument('-d', '--del_cost', default=1, type=float, metavar='',
                        help='Cost of deletion. Must be seen in relation to the specified substitution matrix.')
    parser.add_argument('-c', '--ins_cost', default=1, type=float, metavar='',
                        help='Cost of insertion. Must be seen in relation to the specified substitution matrix.')
    parser.add_argument('-z', '--fig_size', default='20,20,10', type=str, metavar='',
                        help="Figure size string, encoded as 'x,y, labels_font_size'.")

    args = parser.parse_args()
    if args.v_gene is not None and args.v_gene_cost is None:
        raise ValueError('Both, or none, of --v_gene and --v_gene_cost must be specified')
    if args.j_gene is not None and args.j_gene_cost is None:
        raise ValueError('Both, or none, of --j_gene and --j_gene_cost must be specified')
    return args

def possibly_show_matrices_for_argument(args, logger):
    if args.infile == 'show' and args.sequence_col == 'matrices':
        logger.log('identity_matrix')
        logger.log('fast_levenshtein (basic levenshtein, Cython implementation)')
        for m in biomat.available_matrices:
            logger.log(str(m))
        logger.log('peptide_mhc_matrix (from  doi: 10.1186/1471-2105-10-394)')
        exit(0)

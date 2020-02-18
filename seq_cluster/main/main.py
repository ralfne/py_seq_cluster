from seq_cluster.main.arguments import create_arguments, possibly_show_matrices_for_argument
from Logger import StdOutLogger

#https://stats.stackexchange.com/questions/398336/fixing-the-maximum-distance-within-a-cluster

#C:/tmp/seq_cluster/iris_rep.txt "Chain: TRB (1)" -v "TRB - V gene (1)" -j "TRB - J gene (1)" -y 15 -i 10 -n "Clonotype ID" -n "Donor (Name)"
#C:/tmp/seq_cluster/example.txt Seq -n ID
#C:/tmp/seq_cluster/iris_rep.txt "Chain: TRB (1)" -v "TRB - V gene (1)" -j "TRB - J gene (1)" -y 2 -i 1 -n "Clonotype ID" -n "Donor (Name)" -c 1.1 -d 1.1 -m blosum62

#C:/tmp/seq_cluster/iris_rep_subset.txt "Chain: TRB (1)" -v "TRB - V gene (1)" -j "TRB - J gene (1)" -y 15 -i 10 -n "Clonotype ID" -n "Donor (Name)" -n "TRB - V gene (1)" -n "TRB - J gene (1)" -z "40,40,4"
#C:/tmp/seq_cluster/iris_rep.txt "Chain: TRB (1)" -v "TRB - V gene (1)" -j "TRB - J gene (1)" -y 2 -i 1 -n "Clonotype ID" -n "Donor (Name)" -c 1.1 -d 1.1 -m fast_levenshtein
logger = StdOutLogger(verbose=False)

cmd_args = create_arguments()
if cmd_args.outfile is not None:
    import matplotlib
    matplotlib.use('Agg')
    logger.log('matplotlib: using "Agg"...', includeTimestamp=True, onlyIfVerbose=False)

possibly_show_matrices_for_argument(cmd_args, logger)
logger.log('Starting sequence clustering...', includeTimestamp=True, onlyIfVerbose=False)
from seq_cluster.main.execution import run
run(cmd_args, logger)
logger.log('Sequence clustering done.', includeTimestamp=True,onlyIfVerbose=False)

from seq_cluster.clustering.pairwise_distance_matrix import PairwiseDistanceMatrix
from seq_cluster.clustering.simple_abs_distance_calculator import SimpleAbsDistanceCalculator
from seq_cluster.clustering.sequence_distance_calculator import SequenceDistanceCalculator
from seq_cluster.data.data_set import Dataset
from seq_cluster.visualization.clustermap_visualizer import ClustermapVisualizer
import matplotlib.pyplot as plt


fn = 'C:/tmp/seq_cluster/example.txt'

ds = Dataset(fn, sequence_col='Seq', name_cols=['ID', 'Name'], aux_col='Size')

seq_pairs, aux_pairs = ds.get_pairwise()

seq_calc=SequenceDistanceCalculator()
aux_calc = SimpleAbsDistanceCalculator()

seq_matrix = PairwiseDistanceMatrix(seq_pairs, seq_calc)
aux_matrix = PairwiseDistanceMatrix(aux_pairs, aux_calc)

seq_df = seq_matrix.get_as_dataframe()
aux_df = aux_matrix.get_as_dataframe()

viz = ClustermapVisualizer(seq_df, 'average', aux_df,'average')
#viz = ClustermapVisualizer(seq_df, 'average')
viz.run(None)#'C:/tmp/seq_cluster/fig.png')
plt.show()
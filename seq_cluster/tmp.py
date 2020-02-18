from Logger import StdOutLogger

from seq_cluster.distances.bio_matrix import BioMatrix
from seq_cluster.distances.fully_weighed_levenshtein import FullyWeightedLevenshtein, CharacterInsDel
from seq_cluster.distances.identity_distance_matrix import IdentityDistanceMatrix
from seq_cluster.distances.peptide_mhc_matrix import PeptideMhcMatrix

insdel = CharacterInsDel(1, 1)
subst_matrix = BioMatrix('blosum62')
logger=StdOutLogger()
fwl = FullyWeightedLevenshtein(subst_matrix, insdel,logger=logger)
logger.log('start', includeTimestamp=True)
for i in range(100000):
    x1 = fwl.distance('ASSRTDT', 'ASSRDWQ')
    x2 = fwl.distance('RSRTDT', 'ASSYYYWRDWQ')

logger.log('end', includeTimestamp=True)
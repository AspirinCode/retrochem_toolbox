# The 'lib' folder includes libraries from the following repositories:
# --------------------------------------------------------------------
#   1. SA_Score: https://github.com/rdkit/rdkit/tree/master/Contrib/SA_Score
#   2. SC_Score: https://github.com/connorcoley/scscore

from .SA_Score import sascorer
from .SC_Score.scscore.standalone_model_numpy import SCScorer

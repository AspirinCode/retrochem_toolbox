from rdkit.Chem import AllChem
from rdkit.Chem.BRICS import *
from chemical_compounds.analysis_utils import CompoundScoreUtils


if __name__ == '__main__':
    m = AllChem.MolFromSmiles('CCCOCC')
    m2 = BreakBRICSBonds(m, [((1, ), ('1', ))])
    print(AllChem.MolToSmiles(m2))
    print(CompoundScoreUtils.calculate_sc_score("CCOBr", sc_model="2048_bool", verbose=True))

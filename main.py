from chemical_compounds.analysis_utils import ScoreUtils


if __name__ == '__main__':
    print(ScoreUtils.calculate_sc_score("CCOBr", sc_model="2048_bool", verbose=True))

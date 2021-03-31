from chemical_compounds import CompoundRepresentationConversionUtils
from chemical_compounds import MolecularFingerprintUtils

if __name__ == '__main__':

    mol_smiles = "CCOBr"
    mol = CompoundRepresentationConversionUtils.string_to_mol(mol_smiles)
    ecfp = MolecularFingerprintUtils.construct_ecfp(mol, radius=2, bits=32, output_type="np_array")
    hsfp = MolecularFingerprintUtils.construct_hsfp(mol, radius=2, bits=32, from_atoms=[0, 1])

    print(mol_smiles)
    print(mol)
    print(ecfp)
    print(hsfp)

    #m = AllChem.MolFromSmiles('CCCOCC')
    #m2 = BreakBRICSBonds(m, [((1, ), ('1', ))])
    #print(AllChem.MolToSmiles(m2))
    #print(CompoundScoreUtils.calculate_sc_score("CCOBr", sc_model="2048_bool", verbose=True))

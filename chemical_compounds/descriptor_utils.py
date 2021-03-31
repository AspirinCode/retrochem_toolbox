""" Author: Haris Hasic, Phd Student @ Ishida Lab, Department of Computer Science, Tokyo Institute of Technology """

import numpy as np

from typing import List, Set, Tuple, Union

from rdkit.Chem.AllChem import Mol
from rdkit.Chem.AllChem import DataStructs, GetMorganFingerprintAsBitVect, GetDistanceMatrix

from .analysis_utils import CompoundStructureUtils
from .general_utils import CompoundRepresentationConversionUtils


class MolecularFingerprintUtils:
    """ Description: Group of methods for the handling of chemical compound fingerprint descriptors. """

    @staticmethod
    def np_array_to_binary_vector(np_array: np.ndarray) -> DataStructs.ExplicitBitVect:
        """ Description: Convert a NumPy array to a RDKit ExplicitBitVector. """

        binary_vector = DataStructs.ExplicitBitVect(len(np_array))
        binary_vector.SetBitsFromList(np.where(np_array)[0].tolist())

        return binary_vector

    @staticmethod
    def construct_ecfp(compound: Union[str, Mol], radius: int, bits: int, from_atoms=None, output_type="bit_vector",
                       as_type="np_int") -> Union[DataStructs.ExplicitBitVect, np.ndarray]:
        """ Description: Return the Extended Connectivity Fingerprint (ECFP) representation of the whole chemical
                         compound (default) or just from specific atoms, if it is specified through the 'from_atoms'
                         parameter. The type of the whole fingerprint and the individual bits can be adjusted with the
                         parameters 'output_type' and 'as_type', respectively. """

        if isinstance(compound, str):
            compound = CompoundRepresentationConversionUtils.string_to_mol(compound)

            if compound is None:
                raise Exception("Unable to construct a Mol object from the given SMILES string. "
                                "Please make sure that the SMILES string represents a valid chemical compound.")

        if from_atoms is not None:
            ecfp = GetMorganFingerprintAsBitVect(compound, radius=radius, nBits=bits, fromAtoms=from_atoms)
        else:
            ecfp = GetMorganFingerprintAsBitVect(compound, radius=radius, nBits=bits)

        if output_type == "np_array":
            result_ecfp = np.array([])
            DataStructs.cDataStructs.ConvertToNumpyArray(ecfp, result_ecfp)
            ecfp = result_ecfp.astype(np.int_) if as_type == "np_int" else result_ecfp.astype(np.float_)

        return ecfp

    @staticmethod
    def construct_hsfp(compound: Union[str, Mol], radius: int, bits: int,
                       from_atoms: Union[List[int], Set[int], Tuple[int]], neighbourhood_ext=None) -> np.ndarray:
        """ Description: Return the Hot Spot Fingerprint (HSFP) representation of the whole chemical compound (default)
                         or just from specific atoms, if it is specified through the 'from_atoms' parameter, as a NumPy
                         array. If the parameter 'from_atoms' includes all of the atoms of the molecule, the constructed
                         fingerprint is equivalent to the ECFP representation of the whole molecule. The parameter
                         'neighbourhood_ext' controls how many of the neighbourhood atoms are added to the focus atoms
                         specified by the parameter 'from_atoms'. """

        if isinstance(compound, str):
            compound = CompoundRepresentationConversionUtils.string_to_mol(compound)

            if compound is None:
                raise Exception("Unable to construct a Mol object from the given SMILES string. "
                                "Please make sure that the SMILES string represents a valid chemical compound.")

        distance_matrix = GetDistanceMatrix(compound)

        # Set the weight factor for the generation of the HSFP.
        weight_factor = np.max(distance_matrix) if neighbourhood_ext is None else neighbourhood_ext

        # Generate the base of the HSFP, which is basically the ECFP of the core atoms.
        core_fp = MolecularFingerprintUtils.construct_ecfp(
            compound=compound,
            radius=radius,
            bits=bits,
            from_atoms=from_atoms,
            output_type="np_array"
        )
        hsfp = np.array(core_fp)

        # Iterate through and add other layers to the HSFP.
        for i in range(0, int(weight_factor)):
            # Generate a fingerprint for the distance increment, and calculate the bitwise difference.
            atom_environment = CompoundStructureUtils.get_atom_environment(
                compound=compound,
                mol_atoms=from_atoms,
                n_degree=i + 1
            )
            env_fp = MolecularFingerprintUtils.construct_ecfp(
                compound=compound,
                radius=radius,
                bits=bits,
                from_atoms=atom_environment,
                output_type="np_array"
            )
            diff = np.bitwise_xor(core_fp, env_fp)

            # Add the weighted difference vector to the resulting HSFP vector.
            core_fp = core_fp + diff
            hsfp = hsfp + diff * 1 / (i + 2)

        # Return the resulting HSFP rounded to three decimals.
        return np.round(hsfp, 3)

    @staticmethod
    def tanimoto_similarity(ecfp_a: DataStructs.ExplicitBitVect, ecfp_b: DataStructs.ExplicitBitVect) -> float:
        """ Description: Return the Tanimoto similarity value between two fingerprints. """

        return DataStructs.TanimotoSimilarity(ecfp_a, ecfp_b)

    @staticmethod
    def bulk_tanimoto_similarity(ecfp: DataStructs.ExplicitBitVect,
                                 ecfp_pool: List[DataStructs.ExplicitBitVect]) -> List:
        """ Description: Return the Tanimoto similarity values between a single fingerprint and a pool of
                         fingerprints. """

        return DataStructs.BulkTanimotoSimilarity(ecfp, ecfp_pool)

    @staticmethod
    def dice_similarity(ecfp_a: DataStructs.ExplicitBitVect, ecfp_b: DataStructs.ExplicitBitVect) -> float:
        """ Description: Return the Dice similarity value between two fingerprints. """

        return DataStructs.DiceSimilarity(ecfp_a, ecfp_b)

    @staticmethod
    def bulk_dice_similarity(ecfp: DataStructs.ExplicitBitVect,
                             ecfp_pool: List[DataStructs.ExplicitBitVect]) -> List:
        """ Description: Return the Dice similarity values between a single fingerprint and a pool of fingerprints. """

        return DataStructs.BulkDiceSimilarity(ecfp, ecfp_pool)

    @staticmethod
    def tversky_similarity(ecfp_a: DataStructs.ExplicitBitVect, ecfp_b: DataStructs.ExplicitBitVect,
                           a=0.5, b=1.0) -> float:
        """ Description: Return the Tversky similarity value between two fingerprints using the parameter values
                         a and b. """

        return DataStructs.TverskySimilarity(ecfp_a, ecfp_b, a, b)

    @staticmethod
    def bulk_tversky_similarity(ecfp: DataStructs.ExplicitBitVect, ecfp_pool: List[DataStructs.ExplicitBitVect],
                                a=0.5, b=1.0) -> List:
        """ Description: Return the Tversky similarity values between a single fingerprint and a pool of fingerprints
                         using the parameter values a and b. """

        return DataStructs.BulkTverskySimilarity(ecfp, ecfp_pool, a, b)

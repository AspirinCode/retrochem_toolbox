# Author: Haris Hasic, Phd Student @ Ishida Laboratory, Department of Computer Science, Tokyo Institute of Technology

import os

from typing import List, Set, Tuple, Union

from rdkit.Chem.AllChem import Atom, Bond, Mol
from rdkit.Chem.AllChem import GetDistanceMatrix
from rdkit.Chem.SaltRemover import SaltRemover

from .general_utils import CompoundConversionUtils
from .libs.SA_Score import sascorer
from .libs.SC_Score.standalone_model_numpy import SCScorer


# noinspection PyArgumentList
class CompoundStructureUtils:
    """ Description: Group of methods for the handling of chemical compound structures. """

    @staticmethod
    def get_atom_environment(compound: Union[str, Mol],
                             mol_atoms: Union[List[Atom], Tuple[Atom], List[int], Set[int], Tuple[int]],
                             n_degree=1) -> Set:
        """ Description: Get the indices of all 'n_degree' neighbouring atoms for a collection of Mol atoms or Mol atom
                         indices. The 'compound' parameter can be either a SMILES string or a Mol object. """

        if isinstance(compound, str):
            compound = CompoundConversionUtils.string_to_mol(compound)

            if compound is None:
                raise Exception("Unable to construct a Mol object from the given SMILES string. "
                                "Please make sure that the SMILES string represents a valid chemical compound.")

        mol_atom_indices = []

        for mol_atom in mol_atoms:
            if isinstance(mol_atom, Atom):
                mol_atom_indices.append(mol_atom.GetIdx())
            else:
                mol_atom_indices.append(mol_atom)

        # Input the known atoms in the final result and calculate a distance matrix for the molecule.
        neighbour_indices = [atom_ind for atom_ind in mol_atom_indices]
        distance_matrix = GetDistanceMatrix(compound)

        # Check the distances for all neighbours and add them if they are within the designated distance.
        for atom_ind in mol_atom_indices:
            for ind, dist in enumerate(distance_matrix[atom_ind]):
                if dist <= n_degree:
                    neighbour_indices.append(ind)

        return set(neighbour_indices)

    @staticmethod
    def get_bond_environment(compound: Union[str, Mol],
                             mol_bonds: Union[List[Bond], Tuple[Bond], List[int], Set[int], Tuple[int]],
                             n_degree=1) -> Set:
        """ Description: Get the indices of all 'n_degree' neighbouring atoms for a collection of Mol bonds or Mol bond
                         indices. The 'compound' parameter can be either a SMILES string or a Mol object. """

        if isinstance(compound, str):
            compound = CompoundConversionUtils.string_to_mol(compound)

            if compound is None:
                raise Exception("Unable to construct a Mol object from the given SMILES string. "
                                "Please make sure that the SMILES string represents a valid chemical compound.")

        all_mol_bonds_atom_tuples = []

        for mol_bond in mol_bonds:
            if isinstance(mol_bond, Bond):
                all_mol_bonds_atom_tuples.append((mol_bond.GetBeginAtomIdx(), mol_bond.GetEndAtomIdx()))
            else:
                all_mol_bonds_atom_tuples.append((compound.GetBonds()[mol_bond].GetBeginAtomIdx(),
                                                  compound.GetBonds()[mol_bond].GetEndAtomIdx()))

        all_atom_indices = set([atom_ind for single_bond_atom_tuple in all_mol_bonds_atom_tuples
                                for atom_ind in single_bond_atom_tuple])

        return set(CompoundStructureUtils.get_atom_environment(compound, all_atom_indices, n_degree))

    @staticmethod
    def atom_indices_cover_complete_rings(mol_object: Mol, atom_indices: Union[Set, List, Tuple]) -> bool:
        """ Description: Check if a set of atom indices covers complete aromatic rings. """

        all_rings_indices = mol_object.GetRingInfo().AtomRings()

        for ring_indices in all_rings_indices:
            if set(ring_indices).issubset(set(atom_indices)):
                return True

        return False

    @staticmethod
    def count_atom_index_ring_memberships(mol_object: Mol, atom_index: int) -> int:
        """ Description: Count the number of rings in which the atom with the index 'atom_index' is a member of. """

        all_rings_indices = mol_object.GetRingInfo().AtomRings()

        return len([1 for ring_indices in all_rings_indices if atom_index in ring_indices])

    @staticmethod
    def count_mol_bond_ring_memberships(mol_object: Mol, mol_bond: Union[Bond, int]) -> int:
        """ Description: Count the number of rings in which the bond 'mol_bond' is a member of.
                         The 'mol_bond' parameter can be either a Mol bond or a Mol bond index."""

        all_rings_indices = mol_object.GetRingInfo().AtomRings()

        if isinstance(mol_bond, Bond):
            begin_atom_ind = mol_bond.GetBeginAtomIdx()
            end_atom_ind = mol_bond.GetEndAtomIdx()
        else:
            begin_atom_ind = mol_object.GetBonds()[mol_bond].GetBeginAtomIdx()
            end_atom_ind = mol_object.GetBonds()[mol_bond].GetEndAtomIdx()

        return len([1 for ring_indices in all_rings_indices
                    if begin_atom_ind in ring_indices and end_atom_ind in ring_indices])

    @staticmethod
    def get_rest_of_ring_atoms(mol_object: Mol, atom_indices: Union[Set, List, Tuple]) -> List:
        """ Description: Gets the rest of ring atoms for all aromatic rings atoms covered by 'atom_indices'. """

        all_rings_indices = mol_object.GetRingInfo().AtomRings()
        new_atom_indices = [atom_ind for atom_ind in atom_indices]

        while True:
            # Detect whether some of the atom indices are members of a ring, and add the rest of this ring.
            new_additions = []
            for atom_ind in new_atom_indices:
                for ring_ind in all_rings_indices:
                    if atom_ind in ring_ind and not set(ring_ind).issubset(new_atom_indices):
                        new_additions.extend(ring_ind)

            # If there are no detected rings that are not already in the expanded core, break the loop.
            if len(new_additions) == 0:
                break
            else:
                new_atom_indices.extend(list(set(new_additions)))

        return sorted(list(set(new_atom_indices)))

    @staticmethod
    def remove_salts_from_compound(compound: Union[str, Mol], salts_definition_file_path=None,
                                   verbose=False) -> Union[str, None]:
        """ Description: Remove specified salts from a chemical compound using the RDKit salt stripper. If the
                         'salts_definition_file_path' parameter is not specified, only the default RDKit salts will be
                          removed. The 'compound' parameter can be either a SMILES string or a Mol object. """

        try:
            if isinstance(compound, str):
                compound = CompoundConversionUtils.string_to_mol(compound)

                if compound is None:
                    raise Exception("Unable to construct a Mol object from the given SMILES string. "
                                    "Please make sure that the SMILES string represents a valid chemical compound.")

            salt_remover = SaltRemover(defnFilename=salts_definition_file_path)

            return salt_remover.StripMol(compound)

        except Exception as exc_msg:
            if verbose:
                print("Exception occurred during the stripping of the salts from the specified compound. "
                      "Detailed message:\n{}".format(exc_msg))

            return None


class CompoundScoreUtils:
    """ Description: Group of methods for the handling the score calculation for chemical compound structures.
                     Currently includes: SA_Score, and SC_Score. """

    @staticmethod
    def calculate_sa_score(compound: Union[str, Mol], verbose=False) -> Union[float, None]:
        """ Description: Calculate the SA_Score value for a given chemical compound. """

        if isinstance(compound, str):
            compound = CompoundConversionUtils.string_to_mol(compound)

            if compound is None:
                raise Exception("Unable to construct a Mol object from the given SMILES string. "
                                "Please make sure that the SMILES string represents a valid chemical compound.")

        try:
            return sascorer.calculateScore(compound)

        except Exception as exc_msg:
            if verbose:
                print("Exception occurred during the SA_Score calculation. Detailed message:\n{}".format(exc_msg))

            return None

    @staticmethod
    def calculate_sc_score(compound: Union[str, Mol], sc_model="1024bool", verbose=False) -> Union[float, None]:
        """ Description: Calculate the SC_Score value for a given chemical compound. """

        if isinstance(compound, Mol):
            compound = CompoundConversionUtils.mol_to_string(compound)

            if compound is None:
                raise Exception("Unable to construct a SMILES string from the given Mol object. "
                                "Please make sure that the Mol object is a valid chemical compound.")

        try:
            model = SCScorer()

            if sc_model == "1024_uint8":
                model.restore(os.path.join("chemical_compounds/libs/SC_Score/models", "full_reaxys_model_1024uint8",
                                           "model.ckpt-10654.as_numpy.json.gz"))
            if sc_model == "2048_bool":
                model.restore(os.path.join("chemical_compounds/libs/SC_Score/models", "full_reaxys_model_2048bool",
                                           "model.ckpt-10654.as_numpy.json.gz"), FP_len=2048)
            else:
                model.restore(os.path.join("chemical_compounds/libs/SC_Score/models", "full_reaxys_model_1024bool",
                                           "model.ckpt-10654.as_numpy.json.gz"))

            _, sc_score = model.get_score_from_smi(compound)

            return sc_score

        except Exception as exc_msg:
            if verbose:
                print("Exception occurred during the SC_Score calculation. Detailed message:\n{}".format(exc_msg))

            return None

""" Author: Haris Hasic, Phd Student @ Ishida Lab, Department of Computer Science, Tokyo Institute of Technology """

from collections import defaultdict
from copy import deepcopy
from typing import List, Tuple, Union

from rdkit.Chem.AllChem import Atom, Mol, RWMol

from chemical_reactions.general_utils import ReactionRepresentationConversionUtils
from chemical_compounds.general_utils import CompoundRepresentationConversionUtils
from chemical_compounds.analysis_utils import CompoundStructureUtils
from chemical_compounds.editing_utils import EditableStructureUtils


# noinspection PyArgumentList
class ReactionCoreUtils:
    """ Description: Group of methods for the handling of chemical reaction cores. """

    @staticmethod
    def __comp(node, neigh, visited, vis):
        """ Description: Help merging sub-lists that have common elements. """

        nodes = {node}
        next_node = nodes.pop

        while nodes:
            node = next_node()
            vis(node)
            nodes |= neigh[node] - visited

            yield node

    @staticmethod
    def __merge_common(lists):
        """ Description: Merge all sub-list that have common elements. """

        neigh = defaultdict(set)
        visited = set()

        for each in lists:
            for item in each:
                neigh[item].update(each)

        for node in neigh:
            if node not in visited:
                yield sorted(ReactionCoreUtils.__comp(node, neigh, visited, visited.add))

    @staticmethod
    def __atom_in_core(mol_atom_ind: int, reaction_cores: Union[List, Tuple]) -> bool:
        """ Description: Check if a specific atom is in any of the lists of core atoms. """

        for reaction_core in reaction_cores:
            if mol_atom_ind in reaction_core:
                return True

        return False

    @staticmethod
    def __compound_is_mapped(compound: Union[str, Mol]) -> bool:
        """ Description: Check if a compound contains at least one mapped atom. """

        if isinstance(compound, str):
            return ":" in compound

        else:
            for atom in compound.GetAtoms():
                if atom.GetAtomMapNum() != 0:
                    return True

            return False

    @staticmethod
    def same_neighbourhood_size(compound_a: Mol, mol_atom_a: Union[Atom, int], compound_b: Mol,
                                mol_atom_b: Union[Atom, int]) -> bool:
        """ Description: Check whether the same atoms in two different molecules have the same neighbourhood size. """

        if isinstance(mol_atom_a, int):
            mol_atom_a = compound_a.GetAtomWithIdx(mol_atom_a)
        if isinstance(mol_atom_b, int):
            mol_atom_b = compound_b.GetAtomWithIdx(mol_atom_b)

        if len(mol_atom_a.GetNeighbors()) != len(mol_atom_b.GetNeighbors()):
            return False

        return True

    @staticmethod
    def same_neighbour_atoms(compound_a: Mol, mol_atom_a: Union[Atom, int], compound_b: Mol,
                             mol_atom_b: Union[Atom, int]) -> bool:
        """ Description: Check whether the same atoms in two different molecules have retained the same atoms and atom
                         attributes in their immediate neighbourhood according to reaction mapping numbers. """

        if isinstance(mol_atom_a, int):
            mol_atom_a = compound_a.GetAtomWithIdx(mol_atom_a)
        if isinstance(mol_atom_b, int):
            mol_atom_b = compound_b.GetAtomWithIdx(mol_atom_b)

        neighbourhood_a = [(mol_atom.GetAtomMapNum(), mol_atom.GetSymbol(), mol_atom.GetFormalCharge(),
                            mol_atom.GetNumRadicalElectrons(), mol_atom.GetTotalValence())
                           for mol_atom in mol_atom_a.GetNeighbors()]

        neighbourhood_b = [(mol_atom.GetAtomMapNum(), mol_atom.GetSymbol(), mol_atom.GetFormalCharge(),
                            mol_atom.GetNumRadicalElectrons(), mol_atom.GetTotalValence())
                           for mol_atom in mol_atom_b.GetNeighbors()]

        return sorted(neighbourhood_a) == sorted(neighbourhood_b)

    @staticmethod
    def same_neighbour_bonds(compound_a: Mol, mol_atom_a: Union[Atom, int], compound_b: Mol,
                             mol_atom_b: Union[Atom, int]) -> bool:
        """ Description: Check whether the same atoms in two different molecules have retained the same bonds and bond
                         attributes amongst each other in their immediate neighbourhood. """

        if isinstance(mol_atom_a, int):
            mol_atom_a_ind = mol_atom_a
            mol_atom_a = compound_a.GetAtomWithIdx(mol_atom_a)
        else:
            mol_atom_a_ind = mol_atom_a.GetIdx()

        if isinstance(mol_atom_b, int):
            mol_atom_b_ind = mol_atom_b
            mol_atom_b = compound_b.GetAtomWithIdx(mol_atom_b)
        else:
            mol_atom_b_ind = mol_atom_b.GetIdx()

        neighbourhood_1 = [(atom_ind.GetAtomMapNum(),
                            str(compound_a.GetBondBetweenAtoms(mol_atom_a_ind, atom_ind.GetIdx()).GetBondType()))
                           for atom_ind in mol_atom_a.GetNeighbors()]

        neighbourhood_2 = [(atom_ind.GetAtomMapNum(),
                            str(compound_b.GetBondBetweenAtoms(mol_atom_b_ind, atom_ind.GetIdx()).GetBondType()))
                           for atom_ind in mol_atom_b.GetNeighbors()]

        return sorted(neighbourhood_1) == sorted(neighbourhood_2)

    @staticmethod
    def get_reaction_core_atoms(reaction_smiles: str) -> Tuple[List, List]:
        """ Description: Get the indices of atoms that participate in the reaction for each molecule in the reaction.
                         If the molecule does not contain such atoms, return an empty list. This method is based on the
                         assumption that the mapping is correct and done in a 'complete' fashion. This means that all of
                         the atoms in the reactants are mapped, and the ones that persist in the product have the same
                         mapping number. """

        reactants, _, products = ReactionRepresentationConversionUtils.parse_roles_from_reaction_smiles(
            reaction_smiles=reaction_smiles,
            as_what="mol"
        )

        reactants_core_atoms = [set() for _ in range(len(reactants))]
        products_core_atoms = [set() for _ in range(len(products))]

        for p_ind, product in enumerate(products):
            # Only proceed to investigate products that are atom mapped.
            if ReactionCoreUtils.__compound_is_mapped(product):
                for r_ind, reactant in enumerate(reactants):
                    # Only proceed to investigate reactants that are atom mapped.
                    if ReactionCoreUtils.__compound_is_mapped(reactant):

                        for p_atom in product.GetAtoms():
                            # If there are atoms in the product that are not mapped, add them to the core.
                            if p_atom.GetAtomMapNum() <= 0:
                                products_core_atoms[p_ind].add(p_atom.GetIdx())
                                continue

                            for r_atom in reactant.GetAtoms():
                                # If there are atoms in the reactant that are not mapped, add them to the core.
                                if r_atom.GetAtomMapNum() <= 0:
                                    reactants_core_atoms[r_ind].add(r_atom.GetIdx())
                                    continue

                                # If there are atoms in the reactant and product that have the same atom map number,
                                # but different chemical surroundings, add them to the core.
                                if p_atom.GetAtomMapNum() == r_atom.GetAtomMapNum():
                                    if not ReactionCoreUtils.same_neighbourhood_size(product, p_atom.GetIdx(),
                                                                                     reactant, r_atom.GetIdx()) or \
                                       not ReactionCoreUtils.same_neighbour_atoms(product, p_atom.GetIdx(),
                                                                                  reactant, r_atom.GetIdx()) or \
                                       not ReactionCoreUtils.same_neighbour_bonds(product, p_atom.GetIdx(),
                                                                                  reactant, r_atom.GetIdx()):
                                        reactants_core_atoms[r_ind].add(r_atom.GetIdx())
                                        products_core_atoms[p_ind].add(p_atom.GetIdx())

        return reactants_core_atoms, products_core_atoms

    @staticmethod
    def get_reaction_non_core_atoms(reaction_smiles: str) -> Tuple[List, List]:
        """ Description: Get the atoms of the molecule which are not included in the specified reaction cores. This
                         method is just the inverse of the 'get_reaction_core_atoms' method, and all of the same
                         restrictions apply. """

        reactants, _, products = ReactionRepresentationConversionUtils.parse_roles_from_reaction_smiles(
            reaction_smiles=reaction_smiles,
            as_what="mol"
        )

        reactants_non_core_atoms = [set() for _ in range(len(reactants))]
        products_non_core_atoms = [set() for _ in range(len(products))]

        for p_ind, product in enumerate(products):
            for r_ind, reactant in enumerate(reactants):
                for p_atom in product.GetAtoms():

                    # If there are products that are not mapped, add all of their atoms to the non-core.
                    if not ReactionCoreUtils.__compound_is_mapped(product):
                        products_non_core_atoms[p_ind].add(p_atom.GetIdx())
                        continue

                    for r_atom in reactant.GetAtoms():

                        # If there are reactants that are not mapped, add all of their atoms to the non-core.
                        if not ReactionCoreUtils.__compound_is_mapped(reactant):
                            reactants_non_core_atoms[r_ind].add(r_atom.GetIdx())
                            continue

                        # If there are atoms in the reactant and product that have the same atom map number,
                        # and same chemical surroundings, add them to the core.
                        if p_atom.GetAtomMapNum() == r_atom.GetAtomMapNum():
                            if ReactionCoreUtils.same_neighbourhood_size(product, p_atom.GetIdx(),
                                                                         reactant, r_atom.GetIdx()) and \
                               ReactionCoreUtils.same_neighbour_atoms(product, p_atom.GetIdx(),
                                                                          reactant, r_atom.GetIdx()) and \
                               ReactionCoreUtils.same_neighbour_bonds(product, p_atom.GetIdx(),
                                                                          reactant, r_atom.GetIdx()):
                                reactants_non_core_atoms[r_ind].add(r_atom.GetIdx())
                                products_non_core_atoms[p_ind].add(p_atom.GetIdx())

        return reactants_non_core_atoms, products_non_core_atoms

    @staticmethod
    def get_inverse_atoms(reaction_smiles: str, marked_atoms: Tuple[List, List]) -> Tuple[List, List]:
        """ Description: Return the inverse from the marked atoms for each of the reaction roles. """

        reactants, _, products = ReactionRepresentationConversionUtils.parse_roles_from_reaction_smiles(
            reaction_smiles=reaction_smiles,
            as_what="mol_no_maps"
        )

        reaction_roles = [reactants, products]
        reverse_cores = ([], [])

        for role_ind, reaction_role in enumerate(reaction_roles):
            for mol_ind, mol in enumerate(reaction_role):
                local_reverse = set()
                for atom in mol.GetAtoms():
                    if atom.GetIdx() not in marked_atoms[role_ind][mol_ind]:
                        local_reverse.add(atom.GetIdx())

                reverse_cores[role_ind].append(local_reverse)

        return reverse_cores

    @staticmethod
    def get_connected_index_groups(reaction_smiles: str, reaction_cores: Tuple[List, List]):
        """ Description: Get the list of grouped atom indices. This grouping adds another layer in the standard list of
                         lists format. Ideally, the core indices for a single compound in the reaction should all be
                         connected, but sometimes this is not the case. This function can be called to check for such
                         multi-part cores, and to handle them appropriately. If applied to non-core atoms, this method
                         returns the synthon indices. """

        reactants, _, products = ReactionRepresentationConversionUtils.parse_roles_from_reaction_smiles(
            reaction_smiles=reaction_smiles,
            as_what="mol"
        )

        reaction_roles = [reactants, products]
        role_connections, connected_atoms = [[], []], [[], []]

        # Step 1: Aggregate all of the atoms which are connected to each other.
        for rc_ind, reaction_core in enumerate(reaction_cores):
            for rr_ind, reaction_role in enumerate(reaction_core):
                atom_connections = []
                for ind_a, atom_a in enumerate(reaction_role):
                    for ind_b, atom_b in enumerate(reaction_role):
                        if ind_a != ind_b:
                            if reaction_roles[rc_ind][rr_ind].GetBondBetweenAtoms(atom_a, atom_b) is not None:
                                if [atom_a, atom_b] not in atom_connections and \
                                        [atom_b, atom_a] not in atom_connections:
                                    atom_connections.append([atom_a, atom_b])
                role_connections[rc_ind].append(atom_connections)

        # Step 2: Merge all of the individual connections which share the same indices.
        for rlc_ind, role_connection in enumerate(role_connections):
            [connected_atoms[rlc_ind].append(list(ReactionCoreUtils.__merge_common(rc))) for rc in role_connection]

        final_connected_core_indices_groups = deepcopy(connected_atoms)

        # Step 3: Construct the final final connected core indices groups collection.
        for rc_ind, reaction_core in enumerate(reaction_cores):
            for rr_ind, reaction_role in enumerate(reaction_core):
                for atom in reaction_role:
                    if not ReactionCoreUtils.__atom_in_core(atom, connected_atoms[rc_ind][rr_ind]):
                        final_connected_core_indices_groups[rc_ind][rr_ind].append([atom])

        return final_connected_core_indices_groups

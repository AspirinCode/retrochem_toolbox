""" Author: Haris Hasic, Phd Student @ Ishida Lab, Department of Computer Science, Tokyo Institute of Technology """

from rdkit.Chem.AllChem import RWMol, SanitizeMol


# noinspection PyArgumentList
class EditableStructureUtils:
    """ Description: Group of methods for the handling of editable chemical compound structures. """

    @staticmethod
    def get_atoms_to_remove(editable_mol: RWMol):
        """ Description: Fetch all atoms that were marked for removal by the wildcard symbol '*'. """

        atoms_to_remove = []

        for atom in editable_mol.GetAtoms():
            if atom.GetSymbol() == "*":
                atoms_to_remove.append(atom.GetIdx())

        # Return the descending sorted list of atom indices to avoid errors during the removal of the atoms.
        return sorted(list(set(atoms_to_remove)), reverse=True)

    @staticmethod
    def get_bonds_to_remove(editable_mol: RWMol):
        """ Description: Fetch all bond atoms that were marked for removal by the wildcard symbol '*'. """

        bonds_to_remove = []

        for bond in editable_mol.GetBonds():
            if editable_mol.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetSymbol() == "*" and \
                    editable_mol.GetAtomWithIdx(bond.GetEndAtomIdx()).GetSymbol() == "*":
                bonds_to_remove.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))

        # Return the descending sorted list of bond atom indices to avoid errors during the removal of the bond atoms.
        return sorted(list(set(bonds_to_remove)), key=lambda k: k[0], reverse=True)

    @staticmethod
    def remove_floating_atoms(editable_mol: RWMol):
        """ Description: Remove all wildcard atoms '*' that are disconnected from the rest of the compound atoms. """

        leftover_floating_atoms = sorted([atom.GetIdx() for atom in editable_mol.GetAtoms()
                                          if atom.GetSymbol() == "*" and atom.GetDegree() == 0], reverse=True)

        [editable_mol.RemoveAtom(rm_atom) for rm_atom in leftover_floating_atoms]

    @staticmethod
    def remove_marked_atoms(editable_mol: RWMol):
        """ Description: Remove all atoms that were marked for removal by the wildcard symbol '*'. """

        [editable_mol.RemoveAtom(rm_atom) for rm_atom in EditableStructureUtils.get_atoms_to_remove(editable_mol)]

    @staticmethod
    def remove_marked_bonds(editable_mol: RWMol):
        """ Description: Remove all bonds that were marked for removal by the wildcard symbol '*'. """

        [editable_mol.RemoveBond(rm_bond[0], rm_bond[1])
         for rm_bond in EditableStructureUtils.get_bonds_to_remove(editable_mol)]

        # Clean the editable mol from fully disconnected wildcard atoms '*'.
        EditableStructureUtils.remove_floating_atoms(editable_mol)

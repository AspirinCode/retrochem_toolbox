""" Author: Haris Hasic, Phd Student @ Ishida Lab, Department of Computer Science, Tokyo Institute of Technology """

from collections import defaultdict
from copy import deepcopy
from typing import List, Tuple, Union

from rdkit.Chem.AllChem import Atom, Mol, RWMol, SanitizeMol

from chemical_reactions.general_utils import ReactionRepresentationConversionUtils
from chemical_reactions.analysis_utils import ReactionCoreUtils
from chemical_compounds.general_utils import CompoundRepresentationConversionUtils
from chemical_compounds.analysis_utils import CompoundStructureUtils
from chemical_compounds.editing_utils import EditableStructureUtils


# noinspection PyArgumentList
class ReactionAnalysisUtils:
    """ Description: Group of methods for the handling analysis of chemical reaction information. """

    # noinspection PyBroadException
    @staticmethod
    def extract_reactive_substructures(compound: Union[str, Mol], reactive_atoms: List[int]) -> Tuple:
        """ Description: Mark and remove non-reactive atoms from a chemical compound that is participating in a
                         chemical reaction. """

        if isinstance(compound, str):
            compound = CompoundRepresentationConversionUtils.string_to_mol(compound)

            if compound is None:
                raise Exception("Unable to construct a Mol object from the given SMILES string. "
                                "Please make sure that the SMILES string represents a valid chemical compound.")

        editable_mol = RWMol(compound)
        default_editable_mol = deepcopy(editable_mol)

        if len(reactive_atoms) == 0 or len(reactive_atoms) == len(compound.GetAtoms()):
            return editable_mol, default_editable_mol

        nr_atoms = sorted([atom.GetIdx() for atom in compound.GetAtoms() if atom.GetIdx() not in reactive_atoms],
                          reverse=True)

        # Try to remove only atoms that are marked as non-reactive.
        try:
            for rm_atom in nr_atoms:
                editable_mol.ReplaceAtom(rm_atom, Atom("*"))

            editable_mol_ckpt = deepcopy(editable_mol)
            EditableStructureUtils.remove_marked_bonds(editable_mol)
            SanitizeMol(editable_mol)

        # If this fails, usually due to sanitization errors, remove molecule atoms according to the expanded core. Some
        # cores are complete aromatic rings or fused rings which cannot be broken and converted to a RDKit Mol object.
        # Expand the core indices to include all rings or fused rings, and generate a list of atoms that are not part
        # of this core.
        except:
            expanded_core = CompoundStructureUtils.get_rest_of_ring_atoms(compound, reactive_atoms)
            expanded_synthon_atoms = sorted([atom.GetIdx() for atom in compound.GetAtoms()
                                             if atom.GetIdx() not in expanded_core], reverse=True)

            # Create a new copy of the molecule because the previous one may have been modified.
            editable_mol = deepcopy(default_editable_mol)

            # To prevent later sanitization errors due to incorrect valences, skip any isolated atoms that are
            # connected to the core.
            skip_atoms = []
            for bond in editable_mol.GetBonds():
                if bond.GetBeginAtomIdx() in expanded_core and bond.GetEndAtomIdx() not in expanded_core and \
                        editable_mol.GetAtomWithIdx(bond.GetEndAtomIdx()).GetDegree() == 1:
                    skip_atoms.append(bond.GetEndAtomIdx())
                if bond.GetEndAtomIdx() in expanded_core and bond.GetBeginAtomIdx() not in expanded_core and \
                        editable_mol.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetDegree() == 1:
                    skip_atoms.append(bond.GetBeginAtomIdx())

            for atom in set(expanded_synthon_atoms).difference(set(skip_atoms)):
                editable_mol.ReplaceAtom(atom, Atom("*"))

            editable_mol_ckpt = deepcopy(editable_mol)
            EditableStructureUtils.remove_marked_bonds(editable_mol)
            SanitizeMol(editable_mol)

        return editable_mol, editable_mol_ckpt

    # noinspection PyBroadException
    @staticmethod
    def extract_non_reactive_reactant_substructures(reactant_compound: Union[str, Mol],
                                                    reactive_atoms: List[int]) -> Tuple:
        """ Description: Mark and remove reactive atoms from a reactant compound that is participating in a
                         chemical reaction. """

        if isinstance(reactant_compound, str):
            reactant_compound = CompoundRepresentationConversionUtils.string_to_mol(reactant_compound)

            if reactant_compound is None:
                raise Exception("Unable to construct a Mol object from the given SMILES string. "
                                "Please make sure that the SMILES string represents a valid chemical compound.")

        editable_mol = RWMol(reactant_compound)
        default_editable_mol = deepcopy(editable_mol)

        if len(reactive_atoms) == 0 or len(reactive_atoms) == len(reactant_compound.GetAtoms()):
            return editable_mol, default_editable_mol

        # First, try to just remove all of the core atoms without any further modifications.
        try:
            for rm_atom in reactive_atoms:
                editable_mol.ReplaceAtom(rm_atom, Atom("*"))

            editable_mol_ckpt = deepcopy(editable_mol)
            EditableStructureUtils.remove_marked_bonds(editable_mol)

        # If that fails, it's most likely due to the incorrect decomposition of aromatic rings.
        except:
            # Create a new copy of the molecule because the previous one may have been modified.
            editable_mol = deepcopy(default_editable_mol)

            # Try removing all atoms in the present aromatic ring, except atoms shared between fused rings.
            try:
                for bond in editable_mol.GetBonds():
                    if bond.GetBeginAtomIdx() in reactive_atoms and bond.GetEndAtomIdx() in reactive_atoms:
                        if CompoundStructureUtils.count_atom_index_ring_memberships(
                                editable_mol,
                                bond.GetBeginAtomIdx()
                        ) < 2:
                            editable_mol.ReplaceAtom(bond.GetBeginAtomIdx(), Atom("*"))

                        if CompoundStructureUtils.count_atom_index_ring_memberships(
                                editable_mol,
                                bond.GetEndAtomIdx()
                        ) < 2:
                            editable_mol.ReplaceAtom(bond.GetEndAtomIdx(), Atom("*"))

                editable_mol_ckpt = deepcopy(editable_mol)
                EditableStructureUtils.remove_marked_bonds(editable_mol)

            # If this also fails, only remove non-aromatic bond atoms attached to the ring, if there are any.
            except:
                # Create a new copy of the molecule because the previous one may have been modified.
                editable_mol = deepcopy(default_editable_mol)

                for bond in editable_mol.GetBonds():
                    if str(bond.GetBondType()) != "AROMATIC":
                        if bond.GetBeginAtomIdx() in reactive_atoms and bond.GetEndAtomIdx() not in reactive_atoms:
                            editable_mol.ReplaceAtom(bond.GetBeginAtomIdx(), Atom("*"))
                            editable_mol.ReplaceAtom(bond.GetEndAtomIdx(), Atom("*"))
                        elif bond.GetEndAtomIdx() in reactive_atoms and bond.GetBeginAtomIdx() not in reactive_atoms:
                            editable_mol.ReplaceAtom(bond.GetBeginAtomIdx(), Atom("*"))
                            editable_mol.ReplaceAtom(bond.GetEndAtomIdx(), Atom("*"))

                editable_mol_ckpt = deepcopy(editable_mol)
                EditableStructureUtils.remove_marked_bonds(editable_mol)

        return editable_mol, editable_mol_ckpt

    # noinspection PyBroadException
    @staticmethod
    def extract_non_reactive_product_substructures(product_compound: Union[str, Mol],
                                                   reactive_atoms: List[int]) -> Tuple:
        """ Description: Mark and remove reactive atoms from a product compound that is participating in a
                         chemical reaction. """

        if isinstance(product_compound, str):
            product_compound = CompoundRepresentationConversionUtils.string_to_mol(product_compound)

            if product_compound is None:
                raise Exception("Unable to construct a Mol object from the given SMILES string. "
                                "Please make sure that the SMILES string represents a valid chemical compound.")

        editable_mol = RWMol(product_compound)
        default_editable_mol = deepcopy(editable_mol)

        if len(reactive_atoms) == 0 or len(reactive_atoms) == len(product_compound.GetAtoms()):
            return editable_mol

        # First, check if all of the core bonds are aromatic.
        if all([str(bond.GetBondType()) == "AROMATIC" for bond in editable_mol.GetBonds()
                if bond.GetBeginAtomIdx() in reactive_atoms and bond.GetEndAtomIdx() in reactive_atoms]):
            # ------------------------------------------------
            # Reactive atoms CONTAIN ONLY FULL AROMATIC RINGS.
            # ------------------------------------------------
            if CompoundStructureUtils.atom_indices_cover_complete_rings(editable_mol, reactive_atoms):
                # Try removing all aromatic atoms in the core ring, except atoms shared between fused rings.
                try:
                    for bond in editable_mol.GetBonds():
                        if bond.GetBeginAtomIdx() in reactive_atoms and bond.GetEndAtomIdx() in reactive_atoms:
                            if CompoundStructureUtils.count_atom_index_ring_memberships(
                                    editable_mol,
                                    bond.GetBeginAtomIdx()
                            ) < 2:
                                editable_mol.ReplaceAtom(bond.GetBeginAtomIdx(), Atom("*"))

                            if CompoundStructureUtils.count_atom_index_ring_memberships(
                                    editable_mol,
                                    bond.GetEndAtomIdx()
                            ) < 2:
                                editable_mol.ReplaceAtom(bond.GetEndAtomIdx(), Atom("*"))

                    EditableStructureUtils.remove_marked_bonds(editable_mol)

                # If this fails, it is due to kekulization issues. Remove all non-aromatic atoms attached to the ring.
                except:
                    try:
                        # Create a new copy of the molecule because the previous one may have been modified.
                        editable_mol = deepcopy(default_editable_mol)

                        for bond in editable_mol.GetBonds():
                            if str(bond.GetBondType()) != "AROMATIC":
                                if bond.GetBeginAtomIdx() in reactive_atoms and \
                                        bond.GetEndAtomIdx() not in reactive_atoms:
                                    editable_mol.ReplaceAtom(bond.GetBeginAtomIdx(), Atom("*"))
                                    editable_mol.ReplaceAtom(bond.GetEndAtomIdx(), Atom("*"))
                                elif bond.GetEndAtomIdx() in reactive_atoms and \
                                        bond.GetBeginAtomIdx() not in reactive_atoms:
                                    editable_mol.ReplaceAtom(bond.GetBeginAtomIdx(), Atom("*"))
                                    editable_mol.ReplaceAtom(bond.GetEndAtomIdx(), Atom("*"))

                        EditableStructureUtils.remove_marked_bonds(editable_mol)

                    # If that also fails, completely remove all aromatic atoms from the ring.
                    except:
                        # Create a new copy of the molecule because the previous one may have been modified.
                        editable_mol = deepcopy(default_editable_mol)

                        for bond in editable_mol.GetBonds():
                            if str(bond.GetBondType()) == "AROMATIC":
                                if bond.GetBeginAtomIdx() in reactive_atoms and bond.GetEndAtomIdx() in reactive_atoms:
                                    editable_mol.ReplaceAtom(bond.GetBeginAtomIdx(), Atom("*"))
                                    editable_mol.ReplaceAtom(bond.GetEndAtomIdx(), Atom("*"))

                        EditableStructureUtils.remove_marked_bonds(editable_mol)

            # ------------------------------------------------------------------
            # Reactive atoms CONTAIN ALL AROMATIC BONDS, BUT NOT COMPLETE RINGS.
            # ------------------------------------------------------------------
            else:
                # Try removing all of the atoms normally.
                try:
                    for rm_atom in reactive_atoms:
                        editable_mol.ReplaceAtom(rm_atom, Atom("*"))

                    EditableStructureUtils.remove_marked_bonds(editable_mol)

                # If this fails, mark and remove only the nearest non-aromatic bond connected to the aromatic part.
                except:
                    # Create a new copy of the molecule because the previous one may have been modified.
                    editable_mol = deepcopy(default_editable_mol)

                    for bond in editable_mol.GetBonds():
                        if bond.GetBeginAtomIdx() in reactive_atoms and bond.GetEndAtomIdx() in reactive_atoms:
                            if str(bond.GetBondType()) != "AROMATIC":
                                editable_mol.ReplaceAtom(bond.GetBeginAtomIdx(), Atom("*"))
                                editable_mol.ReplaceAtom(bond.GetEndAtomIdx(), Atom("*"))

                    EditableStructureUtils.remove_marked_bonds(editable_mol)

        # -----------------------------------------
        # Reactive atoms CONTAIN NO AROMATIC BONDS.
        # -----------------------------------------
        elif not any([str(bond.GetBondType()) == "AROMATIC" for bond in editable_mol.GetBonds()
                      if bond.GetBeginAtomIdx() in reactive_atoms and bond.GetEndAtomIdx() in reactive_atoms]):
            # Try removing all of the atoms normally. This should work for all cases.
            for rm_atom in reactive_atoms:
                editable_mol.ReplaceAtom(rm_atom, Atom("*"))

            EditableStructureUtils.remove_marked_bonds(editable_mol)

        # -------------------------------------------------------------
        # Reactive atoms CONTAIN COMPLETE RINGS AND NON-AROMATIC BONDS.
        # -------------------------------------------------------------
        else:
            if CompoundStructureUtils.atom_indices_cover_complete_rings(editable_mol, reactive_atoms):
                # Try removing all aromatic atoms in the core ring, except atoms shared between fused rings.
                try:
                    for bond in editable_mol.GetBonds():
                        if bond.GetBeginAtomIdx() in reactive_atoms and bond.GetEndAtomIdx() in reactive_atoms:
                            if CompoundStructureUtils.count_atom_index_ring_memberships(
                                    editable_mol,
                                    bond.GetBeginAtomIdx()
                            ) < 2:
                                editable_mol.ReplaceAtom(bond.GetBeginAtomIdx(), Atom("*"))

                            if CompoundStructureUtils.count_atom_index_ring_memberships(
                                    editable_mol,
                                    bond.GetEndAtomIdx()
                            ) < 2:
                                editable_mol.ReplaceAtom(bond.GetEndAtomIdx(), Atom("*"))

                    EditableStructureUtils.remove_marked_bonds(editable_mol)

                # If this fails, it is  due to kekulization issues. Remove all non-aromatic atoms attached to the ring.
                except:
                    # Create a new copy of the molecule because the previous one may have been modified.
                    editable_mol = deepcopy(default_editable_mol)

                    for bond in editable_mol.GetBonds():
                        if bond.GetBeginAtomIdx() in reactive_atoms and bond.GetEndAtomIdx() in reactive_atoms:
                            if str(bond.GetBondType()) != "AROMATIC":
                                editable_mol.ReplaceAtom(bond.GetBeginAtomIdx(), Atom("*"))
                                editable_mol.ReplaceAtom(bond.GetEndAtomIdx(), Atom("*"))

                    # Sanitize the molecule at this point because some of these entries require it here.
                    SanitizeMol(editable_mol)
                    EditableStructureUtils.remove_marked_bonds(editable_mol)

            # ------------------------------------------------------------------------------
            # Reactive atoms CONTAIN AROMATIC AND NON-AROMATIC BONDS, BUT NO COMPLETE RINGS.
            # ------------------------------------------------------------------------------
            else:
                try:
                    # Try removing all of the atoms normally. This should work for all cases.
                    for rm_atom in reactive_atoms:
                        editable_mol.ReplaceAtom(rm_atom, Atom("*"))

                    EditableStructureUtils.remove_marked_bonds(editable_mol)

                # If this fails, remove only the non-aromatic bonds.
                except:
                    # Create a new copy of the molecule because the previous one may have been modified.
                    editable_mol = deepcopy(default_editable_mol)

                    for bond in editable_mol.GetBonds():
                        if bond.GetBeginAtomIdx() in reactive_atoms and bond.GetEndAtomIdx() in reactive_atoms:
                            if str(bond.GetBondType()) != "AROMATIC":
                                editable_mol.ReplaceAtom(bond.GetBeginAtomIdx(), Atom("*"))
                                editable_mol.ReplaceAtom(bond.GetEndAtomIdx(), Atom("*"))

                    # Remove all bonds that were replaced with the wildcard atom '*'.
                    # If this also fails because of wrongly marked AROMATIC, return default molecule.
                    try:
                        EditableStructureUtils.remove_marked_bonds(editable_mol)
                    except:
                        editable_mol = deepcopy(default_editable_mol)

        return editable_mol, default_editable_mol

    # noinspection PyBroadException
    @staticmethod
    def generate_fragment_data(editable_mol: RWMol, reaction_side="product", basic_editable_mol=None):
        """ Description: Generate various formats of the fragmented molecules. """

        focus_mol = deepcopy(editable_mol)
        smiles = CompoundRepresentationConversionUtils.mol_to_string(focus_mol)
        mol_from_smiles = CompoundRepresentationConversionUtils.string_to_mol(smiles)

        if reaction_side == "reactant":
            # If the editable reactant molecule is split into multiple parts, fall back to the default molecule because
            # for the reactant molecules, just the reactive part needs to be marked.
            if "." in smiles:
                focus_mol = deepcopy(basic_editable_mol)
                smiles = CompoundRepresentationConversionUtils.mol_to_string(focus_mol)
                mol_from_smiles = CompoundRepresentationConversionUtils.string_to_mol(smiles)

            smarts = CompoundRepresentationConversionUtils.mol_to_string(focus_mol, str_format="smarts")
            mol_from_smarts = CompoundRepresentationConversionUtils.string_to_mol(smiles, str_format="smarts")

            if mol_from_smarts is None:
                mol_from_smarts = deepcopy(mol_from_smiles)

            return smiles, smarts, mol_from_smiles, mol_from_smarts

        elif reaction_side == "product":
            all_frag_smiles, all_frag_smarts, all_frag_smiles_mols, all_frag_smarts_mols = [], [], [], []

            for frag_smi in smiles.split("."):
                focus_mol = CompoundRepresentationConversionUtils.string_to_mol(frag_smi)

                all_frag_smiles_mols.append(focus_mol)
                all_frag_smiles.append(CompoundRepresentationConversionUtils.mol_to_string(focus_mol))
                all_frag_smarts.append(CompoundRepresentationConversionUtils.mol_to_string(focus_mol,
                                                                                           str_format="smarts"))

                mol_from_smarts = CompoundRepresentationConversionUtils.string_to_mol(frag_smi, str_format="smarts")

                if mol_from_smarts is None:
                    mol_from_smarts = deepcopy(focus_mol)
                    all_frag_smarts_mols.append(mol_from_smarts)
                    continue
                else:
                    all_frag_smarts_mols.append(mol_from_smarts)

            # Sort the fragments based on the number of atoms in the molecule.
            all_frag_smiles_mols, all_frag_smiles, all_frag_smarts, all_frag_smarts_mols = \
                zip(*sorted(zip(all_frag_smiles_mols, all_frag_smiles, all_frag_smarts, all_frag_smarts_mols),
                            key=lambda k: len(k[0].GetAtoms()), reverse=True))

            # Return the generated data in various formats.
            return list(all_frag_smiles), list(all_frag_smarts), list(all_frag_smiles_mols), list(all_frag_smarts_mols)

        else:
            raise Exception("The only acceptable keywords are 'reactant' and 'product'.")

    @staticmethod
    def extract_info_from_molecule(compound: Union[str, Mol], reactive_atoms: List[int], role="product"):
        """ Description: TBD. """

        if isinstance(compound, str):
            compound = CompoundRepresentationConversionUtils.string_to_mol(compound)

            if compound is None:
                raise Exception("Unable to construct a Mol object from the given SMILES string. "
                                "Please make sure that the SMILES string represents a valid chemical compound.")

        reactive_atoms = sorted(reactive_atoms, reverse=True)

        if role == "reactant":
            rw_mol, basic_rw_mol = ReactionAnalysisUtils.extract_reactive_substructures(compound, reactive_atoms)
            reactive_part = ReactionAnalysisUtils.generate_fragment_data(rw_mol,
                                                                         reaction_side="reactant",
                                                                         basic_editable_mol=basic_rw_mol)

            rw_mol, basic_rw_mol = ReactionAnalysisUtils.extract_non_reactive_reactant_substructures(compound,
                                                                                                     reactive_atoms)
            non_reactive_part = ReactionAnalysisUtils.generate_fragment_data(rw_mol,
                                                                             reaction_side="reactant",
                                                                             basic_editable_mol=basic_rw_mol)

            return reactive_part, non_reactive_part

        else:
            rw_mol, _ = ReactionAnalysisUtils.extract_reactive_substructures(compound, reactive_atoms)
            reactive_part = ReactionAnalysisUtils.generate_fragment_data(rw_mol)

            rw_mol, _ = ReactionAnalysisUtils.extract_non_reactive_product_substructures(compound, reactive_atoms)
            non_reactive_part = ReactionAnalysisUtils.generate_fragment_data(rw_mol)

            return reactive_part, non_reactive_part

    @staticmethod
    def extract_info_from_reaction(reaction_smiles: str, reaction_cores=None):
        """ Extract the reactive and non-reactive parts of the reactant and product molecules from the reaction. """

        reactant_fragments, product_fragments = [], []

        # Extract the reactants and products as RDKit Mol objects and find the reaction cores if none are specified.
        reactants, _, products = ReactionRepresentationConversionUtils.parse_roles_from_reaction_smiles(
            reaction_smiles,
            as_what="mol_no_maps")

        if reaction_cores is None:
            reaction_cores = ReactionCoreUtils.get_reaction_core_atoms(reaction_smiles)

        # Extraction of information from the reactant molecules.
        for r_ind, reactant in enumerate(reactants):
            # Sort the core atom indices in descending order to avoid removal conflicts.
            reactive_atoms = sorted(reaction_cores[0][r_ind], reverse=True)

            rw_mol, basic_rw_mol = ReactionAnalysisUtils.extract_reactive_substructures(reactant, reactive_atoms)
            reactive_part = ReactionAnalysisUtils.generate_fragment_data(rw_mol, reaction_side="reactant",
                                                                         basic_editable_mol=basic_rw_mol)

            rw_mol, basic_rw_mol = ReactionAnalysisUtils.extract_non_reactive_reactant_substructures(reactant,
                                                                                                     reactive_atoms)
            non_reactive_part = ReactionAnalysisUtils.generate_fragment_data(rw_mol, reaction_side="reactant",
                                                                             basic_editable_mol=basic_rw_mol)

            reactant_fragments.append((reactive_part, non_reactive_part))

        # Extraction of information from the product molecules.
        for p_ind, product in enumerate(products):
            # Sort the core atom indices in DESC order to avoid removal conflicts.
            reactive_atoms = sorted(reaction_cores[1][p_ind], reverse=True)

            rw_mol, _ = ReactionAnalysisUtils.extract_reactive_substructures(product, reactive_atoms)
            reactive_part = ReactionAnalysisUtils.generate_fragment_data(rw_mol)

            rw_mol, _ = ReactionAnalysisUtils.extract_non_reactive_product_substructures(product, reactive_atoms)
            non_reactive_part = ReactionAnalysisUtils.generate_fragment_data(rw_mol)

            product_fragments.append((reactive_part, non_reactive_part))

        # Return all of the generated data for a single chemical reaction.
        return reactant_fragments, product_fragments

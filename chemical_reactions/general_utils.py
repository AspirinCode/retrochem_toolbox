""" Author: Haris Hasic, Phd Student @ Ishida Lab, Department of Computer Science, Tokyo Institute of Technology """

import re

from typing import Tuple, Union

from rdkit.Chem.AllChem import ChemicalReaction
from rdkit.Chem.AllChem import ReactionFromSmarts, ReactionToSmiles, ReactionToSmarts, SanitizeRxn

from chemical_compounds.general_utils import CompoundRepresentationConversionUtils


class ReactionRepresentationConversionUtils:
    """ Description: Group of methods for the handling of chemical reaction representation conversions. """

    @staticmethod
    def string_to_reaction(reaction_string: str, str_format="smiles", verbose=True) -> Union[ChemicalReaction, None]:
        """ Description: Convert a chemical reaction string representation (e.g. Reaction SMILES) to a ChemicalReaction
                         object. """

        reaction_object = None

        try:
            if str_format in ["smiles", "smarts"]:
                reaction_object = ReactionFromSmarts(reaction_string)
            else:
                raise Exception("Supported string formats are 'smiles' and 'smarts'. Got: '{}'.".format(str_format))

            SanitizeRxn(reaction_object)

            return reaction_object

        except Exception as exc_msg:
            if verbose:
                if reaction_object is None:
                    print("Exception occurred during the conversion of ", end="")
                else:
                    print("Exception occurred during the sanitization of ", end="")

                print("'{}'. Detailed message:\n{}".format(reaction_string, exc_msg))

            return None

    @staticmethod
    def reaction_to_string(reaction_object: ChemicalReaction, str_format="smiles", additional_transform=True,
                           verbose=False) -> Union[str, None]:
        """ Description: Convert a ChemicalReaction object to a string representation. (e.g. Reaction SMILES).The
                         parameter 'additional_transform' represents whether the output will be canonical in the case
                         of Reaction SMILES. """

        try:
            if str_format == "smiles":
                return ReactionToSmiles(reaction_object, canonical=additional_transform)
            elif str_format == "smarts":
                return ReactionToSmarts(reaction_object)
            else:
                raise Exception("Supported string formats are 'smiles' and 'smarts'. Got: '{}'.".format(str_format))

        except Exception as exc_msg:
            if verbose:
                print("Exception occurred during the conversion of the ChemicalReaction object. "
                      "Detailed message: {}".format(exc_msg))

            return None

    @staticmethod
    def string_to_string(reaction_string: str, str_format="smiles", additional_transform=True,
                         verbose=True) -> Union[str, None]:
        """ Description: Convert a chemical reaction string representation (e.g. Reaction SMILES) to a string
                         representation via the conversion to a ChemicalReaction object. The parameter
                         'additional_transform' represents whether the output will be canonical in the case of Reaction
                         SMILES strings. """

        return ReactionRepresentationConversionUtils.reaction_to_string(
            reaction_object=ReactionRepresentationConversionUtils.string_to_reaction(reaction_string=reaction_string,
                                                                                     str_format=str_format,
                                                                                     verbose=verbose),
            str_format=str_format,
            additional_transform=additional_transform,
            verbose=verbose
        )

    @staticmethod
    def parse_roles_from_reaction_smiles(reaction_smiles: str, as_what="smiles") -> Tuple:
        """ Description: Convert a reaction SMILES string to lists of reactant, reagent and product information. """

        # Split the reaction SMILES string by the '>' symbol to obtain the reactants and products. In cases of extended
        # reaction SMILES strings, there can be additional characters on the product side separated by a whitespace. For
        # this reason, the product side string is always additionally split by a whitespace symbol and only the first
        # element is considered to ensure correct parsing for every reaction SMILES variation.
        reactants = reaction_smiles.split(">")[0].split(".")
        agents = reaction_smiles.split(">")[1].split(".")
        products = reaction_smiles.split(">")[2].split(" ")[0].split(".")

        # Return the original reaction role sub-strings including the reaction atom mappings.
        if as_what == "smiles":
            return [r_smiles for r_smiles in reactants if r_smiles != ""], \
                   [a_smiles for a_smiles in agents if a_smiles != ""], \
                   [p_smiles for p_smiles in products if p_smiles != ""]

        # Return the original reaction role sub-strings excluding the reaction atom mappings.
        elif as_what == "smiles_no_maps":
            return [re.sub(r":[-+]?[0-9]+", "", r_smiles) for r_smiles in reactants if r_smiles != ""], \
                   [re.sub(r":[-+]?[0-9]+", "", a_smiles) for a_smiles in agents if a_smiles != ""], \
                   [re.sub(r":[-+]?[0-9]+", "", p_smiles) for p_smiles in products if p_smiles != ""]

        # Return the lists of atom map numbers for each reaction role.
        elif as_what == "atom_maps":
            return [[int(r_atom_maps[1:]) for r_atom_maps in re.findall(r":[-+]?[0-9]+", r_smiles)]
                    for r_smiles in reactants if r_smiles != ""], \
                   [[int(a_atom_maps[1:]) for a_atom_maps in re.findall(r":[-+]?[0-9]+", a_smiles)]
                    for a_smiles in agents if a_smiles != ""], \
                   [[int(p_atom_maps[1:]) for p_atom_maps in re.findall(r":[-+]?[0-9]+", p_smiles)]
                    for p_smiles in products if products != ""]

        # Return the mol object or the canonical version of the SMILES string for each compound in the reaction roles.
        elif as_what in ["mol", "mol_no_maps", "canonical_smiles", "canonical_smiles_no_maps", "smarts",
                         "smarts_no_maps", "canonical_tuple", "canonical_tuple_no_maps"]:
            all_reaction_role_objects = []

            for reaction_role in [reactants, agents, products]:
                single_reaction_role_objects = []

                for rr_smiles in reaction_role:
                    if rr_smiles != "":
                        if as_what in ["mol", "canonical_smiles", "canonical_tuple"]:
                            mol_object = CompoundRepresentationConversionUtils.string_to_mol(rr_smiles)
                        else:
                            mol_object = CompoundRepresentationConversionUtils.string_to_mol(
                                re.sub(r":[-+]?[0-9]+", "", rr_smiles))

                        if as_what in ["mol", "mol_no_maps"]:
                            single_reaction_role_objects.append(mol_object)
                            continue

                        if as_what in ["smarts", "smarts_no_maps"]:
                            single_reaction_role_objects.append(
                                CompoundRepresentationConversionUtils.mol_to_string(mol_object, str_format="smarts"))
                            continue

                        canonical_smiles = CompoundRepresentationConversionUtils.mol_to_string(mol_object)

                        if as_what in ["canonical_smiles", "canonical_smiles_no_maps"]:
                            single_reaction_role_objects.append(canonical_smiles)
                            continue

                        if as_what in ["canonical_tuple", "canonical_tuple_no_maps"]:
                            single_reaction_role_objects.append((canonical_smiles, mol_object))
                            continue

                all_reaction_role_objects.append(single_reaction_role_objects)

            return tuple(all_reaction_role_objects)

        else:
            raise Exception("Unknown parsing type. Select one of the following: 'smiles', 'smiles_no_maps', 'smarts', "
                            "'smarts_no_maps', 'mol', 'mol_no_maps', 'canonical_smiles', 'canonical_smiles_no_maps'.")

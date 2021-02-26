# Author: Haris Hasic, Phd Student @ Ishida Laboratory, Department of Computer Science, Tokyo Institute of Technology

from typing import Union

from rdkit.Chem.AllChem import Mol
from rdkit.Chem.AllChem import MolFromSmiles, MolFromSmarts, MolToSmiles, MolToSmarts, SanitizeMol


class CompoundConversionUtils:
    """ Description: Group of methods for the handling of chemical compound representation conversions. """

    @staticmethod
    def string_to_mol(compound_string: str, str_format="smiles", verbose=False) -> Union[Mol, None]:
        """ Description: Convert a chemical compound string representation to a Mol object. """

        mol_object = None

        try:
            if str_format == "smiles":
                mol_object = MolFromSmiles(compound_string)
            elif str_format == "smarts":
                mol_object = MolFromSmarts(compound_string)
            else:
                raise Exception("Supported string formats are 'smiles' and 'smarts'. Got: '{}'.".format(str_format))

            SanitizeMol(mol_object)

            return mol_object

        except Exception as exc_msg:
            if verbose:
                if mol_object is None:
                    print("Exception occurred during the conversion of ", end="")
                else:
                    print("Exception occurred during the sanitization of ", end="")

                print("'{}'. Detailed message:\n{}".format(compound_string, exc_msg))

            return None

    @staticmethod
    def mol_to_string(mol_object: Mol, str_format="smiles", canonical=True, verbose=False) -> Union[str, None]:
        """ Description: Convert a chemical compound Mol object to a string representation. """

        try:
            if str_format == "smiles":
                return MolToSmiles(mol_object, canonical=canonical)
            elif str_format == "smarts":
                return MolToSmarts(mol_object)
            else:
                raise Exception("Supported string formats are 'smiles' and 'smarts'. Got: '{}'.".format(str_format))

        except Exception as exc_msg:
            if verbose:
                print("Exception occurred during the conversion of the Mol object. "
                      "Detailed message: {}".format(exc_msg))

            return None

    @staticmethod
    def string_to_canonical_string(compound_string: str, in_str_format="smiles", out_str_format="smiles",
                                   verbose=False) -> Union[str, None]:
        """ Description: Canonicalize a chemical compound string representation. """

        return CompoundConversionUtils.mol_to_string(
            CompoundConversionUtils.string_to_mol(compound_string, str_format=in_str_format, verbose=verbose),
            str_format=out_str_format, canonical=True, verbose=verbose)

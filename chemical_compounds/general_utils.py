""" Author: Haris Hasic, Phd Student @ Ishida Lab, Department of Computer Science, Tokyo Institute of Technology """

from typing import Union

from rdkit.Chem.AllChem import Mol
from rdkit.Chem.AllChem import MolFromSmiles, MolFromSmarts, MolToSmiles, MolToSmarts, SanitizeMol


class CompoundRepresentationConversionUtils:
    """ Description: Group of methods for the handling of chemical compound representation conversions. """

    @staticmethod
    def string_to_mol(compound_string: str, str_format="smiles", verbose=True) -> Union[Mol, None]:
        """ Description: Convert a chemical compound string representation (e.g. SMILES) to a Mol object. """

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
    def mol_to_string(mol_object: Mol, str_format="smiles", additional_transform=True,
                      verbose=False) -> Union[str, None]:
        """ Description: Convert a chemical compound Mol object to a string representation (e.g. SMILES). The parameter
                         'additional_transform' represents whether the output will be canonical or isomeric in the case
                         of SMILES and SMARTS strings, respectively. """

        try:
            if str_format == "smiles":
                return MolToSmiles(mol_object, canonical=additional_transform)
            elif str_format == "smarts":
                return MolToSmarts(mol_object, isomericSmiles=additional_transform)
            else:
                raise Exception("Supported string formats are 'smiles' and 'smarts'. Got: '{}'.".format(str_format))

        except Exception as exc_msg:
            if verbose:
                print("Exception occurred during the conversion of the object. Detailed message: {}".format(exc_msg))

            return None

    @staticmethod
    def string_to_string(compound_string: str, in_str_format="smiles", out_str_format="smiles",
                         additional_transform=True, verbose=False) -> Union[str, None]:
        """ Description: Convert a chemical compound string representation (e.g. SMILES) to a string representation via
                         the conversion to a Mol object. The parameter 'additional_transform' represents whether the
                         output will be canonical or isomeric in the case of SMILES and SMARTS strings,
                         respectively. """

        return CompoundRepresentationConversionUtils.mol_to_string(
            mol_object=CompoundRepresentationConversionUtils.string_to_mol(compound_string=compound_string,
                                                                           str_format=in_str_format,
                                                                           verbose=verbose),
            str_format=out_str_format,
            additional_transform=additional_transform,
            verbose=verbose
        )

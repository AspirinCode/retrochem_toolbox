""" Author: Haris Hasic, Phd Student @ Ishida Lab, Department of Computer Science, Tokyo Institute of Technology """

import io

import numpy as np

from cairosvg import svg2png
from PIL import Image
from typing import Dict, List, Tuple, Union

from rdkit.Chem.AllChem import Compute2DCoords, Mol
from rdkit.Chem.Draw import DrawMorganBits
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DSVG, PrepareMolForDrawing

from .general_utils import CompoundRepresentationConversionUtils


class CompoundVisualizationUtils:
    """ Description: Group of methods for the handling of the chemical compound visualizations. """

    @staticmethod
    def __assign_colors_to_indices(indices_subsets: List[List]) -> Tuple[List, Dict]:
        """ Description: Assign different color codes to different subsets of indices. """

        if indices_subsets is None:
            return [], {}
        elif len(indices_subsets) == 0:
            return [], {}

        unified_indices, colors = [], {}

        for indices_subset in indices_subsets:
            indices_subset_color = tuple(np.random.rand(3,))

            for index in indices_subset:
                unified_indices.append(index)
                colors.update({index: indices_subset_color})

        return unified_indices, colors

    # noinspection PyArgumentList
    @staticmethod
    def draw_compound(compound: Union[str, Mol], im_size_x=300, im_size_y=200, highlight_atoms=None,
                      highlight_bonds=None) -> Image:
        """ Description: Draw the chemical compound with or without highlighted individual atoms or bonds. """

        if isinstance(compound, str):
            compound = CompoundRepresentationConversionUtils.string_to_mol(compound)

            if compound is None:
                raise Exception("Unable to construct a Mol object from the given SMILES string. "
                                "Please make sure that the SMILES string represents a valid chemical compound.")

        highlight_atoms, highlight_atom_colors = CompoundVisualizationUtils.__assign_colors_to_indices(highlight_atoms)
        highlight_bonds, highlight_bond_colors = CompoundVisualizationUtils.__assign_colors_to_indices(highlight_bonds)

        Compute2DCoords(compound)

        try:
            compound.GetAtomWithIdx(0).GetExplicitValence()
        except RuntimeError:
            compound.UpdatePropertyCache(False)
        try:
            new_mol = PrepareMolForDrawing(compound, kekulize=True)
        except ValueError:
            new_mol = PrepareMolForDrawing(compound, kekulize=False)

        drawer = MolDraw2DSVG(im_size_x, im_size_y)

        if highlight_atoms is not None and highlight_bonds is None:
            drawer.DrawMolecule(new_mol, highlightAtoms=highlight_atoms, highlightAtomColors=highlight_atom_colors)
        elif highlight_atoms is None and highlight_bonds is not None:
            drawer.DrawMolecule(new_mol, highlightBonds=highlight_bonds, highlightBondColors=highlight_bond_colors)
        else:
            drawer.DrawMolecule(new_mol, highlightAtoms=highlight_atoms, highlightAtomColors=highlight_atom_colors,
                                highlightBonds=highlight_bonds, highlightBondColors=highlight_bond_colors)

        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace("svg:", "")

        return Image.open(io.BytesIO(svg2png(svg)))

    @staticmethod
    def draw_fingerprint_substructures(compound: Union[str, Mol], radius: int, from_atoms=None, im_size_x=300,
                                       im_size_y=300) -> Image:
        """ Description: Draw the fingerprint substructures of a chemical compound for a specified radius. """

        if isinstance(compound, str):
            compound = CompoundRepresentationConversionUtils.string_to_mol(compound)

            if compound is None:
                raise Exception("Unable to construct a Mol object from the given SMILES string. "
                                "Please make sure that the SMILES string represents a valid chemical compound.")

        bit_info = {}
        fp = GetMorganFingerprintAsBitVect(compound, radius=radius, fromAtoms=from_atoms, bitInfo=bit_info)
        on_bits = [(compound, x, bit_info) for x in fp.GetOnBits()]

        drawer = DrawMorganBits(on_bits, molsPerRow=3, subImgSize=(im_size_x, im_size_y),
                                legends=[str(x) for x in fp.GetOnBits()])
        svg = drawer.replace("svg:", "")

        return Image.open(io.BytesIO(svg2png(svg)))

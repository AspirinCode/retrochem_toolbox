""" Author: Haris Hasic, Phd Student @ Ishida Lab, Department of Computer Science, Tokyo Institute of Technology """

from PIL import Image
from typing import Union

from rdkit.Chem.AllChem import ChemicalReaction

from .general_utils import ReactionRepresentationConversionUtils
from chemical_compounds.visualization_utils import CompoundVisualizationUtils


class ReactionVisualizationUtils:
    """ Description: Group of methods for the handling of the chemical compound visualizations. """

    @staticmethod
    def draw_reaction(reaction: Union[str, ChemicalReaction], show_mapping=True, show_reagents=True,
                      reaction_cores=None, im_size_x=300, im_size_y=200) -> Image:
        """ Description: Draw the chemical reaction with or without highlighted reaction cores and reactive parts. """

        if not isinstance(reaction, str):
            reaction = ReactionRepresentationConversionUtils.reaction_to_string(reaction)

            if reaction is None:
                raise Exception("Unable to convert the ChemicalReaction object to a Reaction SMILES string. "
                                "Please make sure that the ChemicalReaction object represents a valid reaction.")

        reactants, reagents, products = ReactionRepresentationConversionUtils.parse_roles_from_reaction_smiles(
            reaction_smiles=reaction,
            as_what="mol" if show_mapping else "mol_no_maps"
        )

        if reaction_cores is None:
            reaction_cores = [[], []]

        mol_images = []

        for r_ind, reactant in enumerate(reactants):
            if len(reaction_cores[0]) > 0:
                mol_images.append(
                    CompoundVisualizationUtils.draw_compound(
                        compound=reactant,
                        im_size_x=im_size_x,
                        im_size_y=im_size_y,
                        highlight_atoms=[reaction_cores[0][r_ind]]
                    )
                )
            else:
                mol_images.append(
                    CompoundVisualizationUtils.draw_compound(
                        compound=reactant,
                        im_size_x=im_size_x,
                        im_size_y=im_size_y
                    )
                )

            if r_ind == len(reactants) - 1:
                mol_images.append(Image.open("chemical_reactions/assets/arrow.png"))
            else:
                mol_images.append(Image.open("chemical_reactions/assets/plus.png"))

        if len(reagents) > 0 and show_reagents:
            for rg_ind, reagent in enumerate(reagents):
                mol_images.append(
                    CompoundVisualizationUtils.draw_compound(
                        compound=reagent,
                        im_size_x=im_size_x,
                        im_size_y=im_size_y
                    )
                )
                if rg_ind == len(reagents) - 1:
                    mol_images.append(Image.open("chemical_reactions/assets/arrow.png"))
                else:
                    mol_images.append(Image.open("chemical_reactions/assets/plus.png"))

        for p_ind, product in enumerate(products):
            if len(reaction_cores[1]) > 0:
                mol_images.append(
                    CompoundVisualizationUtils.draw_compound(
                        compound=product,
                        im_size_x=im_size_x,
                        im_size_y=im_size_y,
                        highlight_atoms=[reaction_cores[1][p_ind]]
                    )
                )
            else:
                mol_images.append(
                    CompoundVisualizationUtils.draw_compound(
                        compound=product,
                        im_size_x=im_size_x,
                        im_size_y=im_size_y
                    )
                )

            if p_ind != len(products) - 1:
                mol_images.append(Image.open("chemical_reactions/assets/plus.png"))

        widths, heights = zip(*(i.size for i in mol_images))
        total_width = sum(widths)
        max_height = max(heights)
        new_im = Image.new("RGB", (total_width, max_height), (255, 255, 255))

        x_offset, y_offset = 0, 0

        for ind, im in enumerate(mol_images):
            if ind % 2 != 0:
                y_offset = round(im_size_y / 2 - im.size[1] / 2)
            else:
                y_offset = 0

            new_im.paste(im, (x_offset, y_offset))
            x_offset += im.size[0]

        return new_im

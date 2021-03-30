import io
from PIL import Image
from cairosvg import svg2png
from rdkit.Chem import AllChem, Draw, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D


class CompoundVisualizationUtils:
    """ Description: Group of methods for the handling of chemical compound structures. """

    @staticmethod
    def __assign_colors_to_indices(indices_subsets):
        """ Assigns different colors to different subsets of indices. """

        # If there are no highlighted elements, return no colors.
        if indices_subsets is None:
            return [], {}

        # Define the colors that will be used for highlighting different groups of elements.
        color_codes = {1: (0.9, 0.4, 0.4), 2: (0.1, 0.9, 0.4), 3: (0.1, 0.4, 0.9), 4: (0.9, 1, 0.4), 5: (0.9, 0.4, 0.9)}
        colors, unified_indices = {}, []

        # Add colors to different subsets.
        color_key = 1
        for subset in indices_subsets:
            for s in subset:
                unified_indices.append(s)

                if color_key in color_codes.keys():
                    colors.update({s: color_codes[color_key]})
                else:
                    colors.update({s: color_codes[color_key - 1]})

            color_key = color_key + 1

        # Return the generated colors.
        return unified_indices, colors
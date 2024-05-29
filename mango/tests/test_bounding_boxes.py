from mango.design_spaces.bounding_box import TriclinicBox, MonoclinicBox, OrthorhombicBox, TetragonalBox, \
    RhombohedralBox, HexagonalBox, CubicBox

from mango.mango_features.preserved_regions import PreservedVertex
from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
from mango.visualizations import display_mesh_design
import numpy as np


if __name__ == '__main__':
    new_box1 = CubicBox(a=50)
    new_box2 = TriclinicBox(a=50, b=60, c=90, alpha=70, beta=80, gamma=90)
    new_box3 = MonoclinicBox(a=50, b=80, c=40, beta=60)
    new_box4 = OrthorhombicBox(a=50, b=20, c=90)
    new_box5 = TetragonalBox(a=50, c=80)
    new_box6 = RhombohedralBox(a=35, alpha=80)
    new_box7 = HexagonalBox(a=25, c=50)

    preserved_regions = []

    excluded_regions = []

    design_space = PolyhedralSpace(bounding_box=new_box1, preserved=preserved_regions, excluded=excluded_regions)

    # Plot and show the design to validate it has been created successfully:
    cylinder_plot = display_mesh_design.CylindricalRepresentation(design=design_space, bounding_box=new_box1)
    #cylinder_plot.show_figure()


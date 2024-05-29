from mango.design_spaces.bounding_box import TetragonalBox
from mango.mango_features.preserved_regions import PreservedVertex
from mango.mango_features.excluded_regions import RectangularPrism,Sphere
from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
from mango.visualizations import display_mesh_design
import numpy as np


if __name__ == '__main__':
    X = 50
    Y = 50
    Z = 70

    new_box = TetragonalBox(a=X, c=Z)
    preserved_regions = [
        PreservedVertex(v1=np.array([X / 2, 0, Z / 2])),
        PreservedVertex(v1=np.array([X / 2, Y, Z / 2])),
        PreservedVertex(v1=np.array([0, Y / 2, Z / 2])),
        PreservedVertex(v1=np.array([X, Y / 2, Z / 2])),
        PreservedVertex(v1=np.array([X / 2, Y / 2, 0])),
        PreservedVertex(v1=np.array([X / 2, Y / 2, Z])),
    ]

    excluded_regions = [
        RectangularPrism(c1=np.array([0, 0, 0]), c2=np.array([15, 15, 15])),
        Sphere(diameter=15, center=np.array([X/2, Y/2, Z/2]))
    ]

    design_space = PolyhedralSpace(bounding_box=new_box, preserved=preserved_regions, excluded=excluded_regions)
    design_space.generate_start_shape()

    # Plot and show the design to validate it has been created successfully:
    mesh_plot = display_mesh_design.MeshRepresentation(design=design_space, bounding_box=new_box)
    mesh_plot.show_figure()
    cylinder_plot = display_mesh_design.CylindricalRepresentation(design=design_space, bounding_box=new_box)
    cylinder_plot.show_figure()


from mango.design_spaces.bounding_box import TetragonalBox
from mango.mango_features.preserved_regions import PreservedEdge
from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
from mango.visualizations import display_mesh_design
import numpy as np


if __name__ == '__main__':
    X = 80
    Y = 80
    Z = 60.

    new_box = TetragonalBox(a=X, c=Z)
    preserved_regions = [
        PreservedEdge(v1=np.array([30., 50., 0]), v2=np.array([30., 30., 0])),
        PreservedEdge(v1=np.array([30., 30., 0]), v2=np.array([50., 30., 0])),
        PreservedEdge(v1=np.array([50., 30., 0]), v2=np.array([50., 50., 0])),
        PreservedEdge(v1=np.array([50., 50., 0]), v2=np.array([30., 50., 0])),
        PreservedEdge(v1=np.array([30., 50., Z]), v2=np.array([30., 30., Z])),
        PreservedEdge(v1=np.array([30., 30., Z]), v2=np.array([50., 30., Z])),
        PreservedEdge(v1=np.array([50., 30., Z]), v2=np.array([50., 50., Z])),
        PreservedEdge(v1=np.array([50., 50., Z]), v2=np.array([30., 50., Z])),
    ]

    excluded_regions = [
    ]

    design_space = PolyhedralSpace(bounding_box=new_box, preserved=preserved_regions, excluded=excluded_regions)
    design_space.generate_start_shape()

    # Plot and show the design to validate it has been created successfully:
    mesh_plot = display_mesh_design.MeshRepresentation(design=design_space, bounding_box=new_box)
    mesh_plot.show_figure()
    cylinder_plot = display_mesh_design.CylindricalRepresentation(design=design_space, bounding_box=new_box,
                                                                  cylinder_diameter=3.75)
    cylinder_plot.show_figure()


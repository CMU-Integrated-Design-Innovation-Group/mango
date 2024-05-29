from mango.design_spaces.bounding_box import CubicBox
from mango.mango_features.preserved_regions import PreservedVertex, PreservedEdge
from mango.mango_features.excluded_regions import Sphere
from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
from mango.visualizations import display_mesh_design
import numpy as np


if __name__ == '__main__':
    X = 50
    new_box = CubicBox(a=X)

    preserved_regions = [
        PreservedVertex(v1=np.array([25, 0, 25])),
        PreservedEdge(v1=np.array([15, 50, 25]), v2=np.array([35, 50, 25])),
        PreservedVertex(v1=np.array([25, 25, 0])),

        PreservedEdge(v1=np.array([0., 15., 0.]), v2=np.array([0., 35., 0.])),

        PreservedEdge(v1=np.array([35., 35., 50]), v2=np.array([35., 15., 50])),
        PreservedEdge(v1=np.array([35., 15., 50]), v2=np.array([15., 15., 50])),
        PreservedEdge(v1=np.array([15., 15., 50]), v2=np.array([15., 35., 50])),
        PreservedEdge(v1=np.array([15., 35., 50]), v2=np.array([35., 35., 50])),
    ]

    excluded_regions = [
        Sphere(diameter=10, center=np.array([12.5, 16.666, 16.666])),
        Sphere(diameter=10, center=np.array([27.5, 37.5, 25])),
    ]

    design_space = PolyhedralSpace(bounding_box=new_box, preserved=preserved_regions, excluded=excluded_regions)
    design_space.generate_start_shape()

    # Plot and show the design to validate it has been created successfully:
    mesh_plot = display_mesh_design.MeshRepresentation(design=design_space, bounding_box=new_box)
    mesh_plot.show_figure()
    cylinder_plot = display_mesh_design.CylindricalRepresentation(design=design_space, bounding_box=new_box)
    cylinder_plot.show_figure()


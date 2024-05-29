"""
A.J. Vetturini
IDIG and MMBL
Carnegie Mellon University

This modules defines preserved region areas that will be used in defining the polyhedral_design_space to ensure DNA
is located in specific places in the design space. A user can specify either an edge of DNA or a vertex which represents
DNA bundles "meeting" at a specific point in the wireframe shape.
"""
# Import Modules
from dataclasses import dataclass
import numpy as np
import networkx as nx

@dataclass
class PreservedVertex(object):
    """
    The preserved vertex will guarantee a specific point in space is held constantly at that location during the
    generative process. Note that the input conditions defining the Preserved Regions MUST be able to be triangulated
    via the alphashape algorithm.

    Parameters
    ------------
    v1 : Numpy array defining where the vertex is located in 3D space
    """
    v1: np.array

    def __post_init__(self):
        self.relabel_dict = {}
        # First we check to see what kind of face we have:
        self.binding_type = 'vertex'
        graph_edges = []

        # Next we create a NetworkX Graph:
        self.preserved_region_as_graph = nx.Graph()  # Initialize a graph that will store data of these points
        self.graph_points = [
            (0, {'x': self.v1[0], 'y': self.v1[1], 'z': self.v1[2], 'terminal': True, 'preserved_region': True})
        ]

        # Now add nodes and points to the NetworkX Graph:
        self.preserved_region_as_graph.add_nodes_from(self.graph_points)
        self.preserved_region_as_graph.add_edges_from(graph_edges)



@dataclass
class PreservedEdge(object):
    """
    The preserved edge will guarantee a bundle of DNA is located connecting two vertices in 3D space. Note that the
    input conditions defining the Preserved Regions MUST be able to be triangulated via the alphashape algorithm.

    Parameters
    ------------
    v1 : Numpy array defining where the 1st vertex is located in 3D space
    v3 : Numpy array defining where the 2nd vertex is located in 3D space
    """
    v1: np.array
    v2: np.array

    def __post_init__(self):
        self.relabel_dict = {}
        # First we check to see what kind of face we have:
        self.binding_type = 'edge'
        v_list = [self.v1, self.v2]
        graph_edges = [(0, 1)]

        # Next we create a NetworkX Graph:
        self.preserved_region_as_graph = nx.Graph()  # Initialize a graph that will store data of these points
        graph_points = []
        for ct, i in enumerate(v_list):
            graph_points.append((ct, {'x': i[0], 'y': i[1], 'z': i[2], 'terminal': True, 'preserved_region': True}))

        # Now add nodes and points to the NetworkX Graph:
        self.preserved_region_as_graph.add_nodes_from(graph_points)
        self.preserved_region_as_graph.add_edges_from(graph_edges)

        # Append a list for all the points found in the preserved_region:
        self.graph_points = []
        for i in graph_points:
            self.graph_points.append(i[1])
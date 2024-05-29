"""
A.J. Vetturini
IDIG and MMBL
Carnegie Mellon University

These are various functions used throughout the application that range from mathematical functions to graph-based
functions which return specific values from a provided graph. There are no custom data types used in any of these
functions which is the criteria I used to "separate" out the code.
"""
import networkx as nx
import numpy as np
from itertools import combinations
from numba import jit

def xyz_from_graph(graph: nx.Graph, node: int) -> np.array:
    """
    This function returns the X Y Z value of a specified node # from the design_graph

    :param graph: NetworkX undirected graph which is stored as the design_graph object
    :param node: Node # to return the X Y Z location of
    """
    pt = np.array([graph.nodes[node]['x'], graph.nodes[node]['y'], graph.nodes[node]['z']])
    return pt


def node_from_xyz(graph: nx.Graph, xyz: np.array) -> int:
    """
    This function returns the node number of a specified X Y Z array from the passed in design graph

    :param graph: NetworkX undirected graph which is stored as the design_graph object
    :param xyz: NumPy array in form of [X, Y, Z] to return the node of
    :return: Integer of node number represented by xyz
    """
    tolerance = 1e-2  # I use a little bit larger tolerance since we usually aren't going to low decimal-wise
                      # due to inherent precision limitations
    for n in graph.nodes:
        #if graph.nodes[n]['x'] == xyz[0] and graph.nodes[n]['y'] == xyz[1] and graph.nodes[n]['z'] == xyz[2]:
        if (np.isclose(graph.nodes[n]['x'], xyz[0], atol=tolerance) and
            np.isclose(graph.nodes[n]['y'], xyz[1], atol=tolerance) and
            np.isclose(graph.nodes[n]['z'], xyz[2], atol=tolerance)):
            return n
    raise Exception('Unable to find node in the design graph!!!')


def nodes_on_same_face(all_faces: dict, n1: int, n2: int) -> bool:
    """
    This function will take in two nodes and return TRUE if they ARE on the same face, and FALSE if not

    :param all_faces: Dictionary of all faces contained in the design
    :param n1: 1st node to check
    :param n2: 2nd node to check
    :return: True: Nodes on same face; False: Nodes NOT on same face.
    """
    for f in all_faces.keys():
        if n1 in list(f) and n2 in list(f):
            # If both nodes are in a tuple that is defining a face, then the edge must exist.
            return True
    # If both nodes are not in any face, then we have a new edge.
    return False


def all_node_xyz_dict(graph: nx.Graph, nodes: list) -> dict:
    """
    This function simply returns all nodes as a dictionary for design constraint checking

    :param graph: NetworkX undirected graph which is stored as the design_graph object
    :param nodes: List of all node index values [1, 2, 3, ...]
    :return: Dictionary mapping node index values to the node X Y Z value
    """
    node_map = {}
    for n in nodes:
        pos = xyz_from_graph(graph=graph, node=n)
        node_map[n] = pos
    return node_map


def length_and_direction_between_nodes(graph: nx.Graph, node1: int, node2: int) -> tuple[np.array, np.array]:
    """
    This function calculates the distance between two nodes and returns that distance and the direction of the vector
    connected N1 -> N2

    :param graph: NetworkX undirected graph which is stored as the design_graph object
    :param node1: Index of node 1 stored in design_graph
    :param node2: Index of node 2 stored in design_graph
    :return: Tuple of (Magnitude / Length between nodes, vector [i, j, k] establishing direction
    """
    P1 = xyz_from_graph(graph=graph, node=node1)
    P2 = xyz_from_graph(graph=graph, node=node2)
    direction = P2 - P1
    return np.linalg.norm(direction), direction


def calculate_design_edge_lengths(graph: nx.Graph) -> list:
    """
    Loop that calculate the total edge lengths of all edges in the design in units of nanometers

    :param graph: NetworkX undirected graph which is stored as the design_graph object
    :return: List of length values [23.23, 24.56, ...] in units (nm)
    """
    all_lengths = []
    for v1, v2 in graph.edges():
        length, _ = length_and_direction_between_nodes(graph=graph, node1=v1, node2=v2)
        all_lengths.append(length)
    return all_lengths


def find_shortest_path(graph: nx.Graph, n1: int, n2: int) -> list:
    """
    This function returns the shortest path between two nodes on the graph

    :param graph: NetworkX undirected graph which is stored as the design_graph object
    :param n1: Index of node 1 stored in design_graph
    :param n2: Index of node 2 stored in design_graph
    :return: List of length values [23.23, 24.56, ...] in units (nm)
    """
    shortest_paths = nx.all_shortest_paths(graph, source=n1, target=n2)
    tot_dists = []
    for shortest_path in shortest_paths:
        tot_dist = 0.
        for i in range(len(shortest_path) - 1):
            current_node = shortest_path[i]
            next_node = shortest_path[i + 1]
            dist, _ = length_and_direction_between_nodes(graph=graph, node1=current_node, node2=next_node)
            tot_dist += dist
        tot_dists.append(tot_dist)
    return tot_dists


def rotation_matrix_from_vectors(vec1: np.array, vec2: np.array) -> np.array:
    """
    Calculates the rotation matrix defined by two vectors.

    If you ever find yourself doing tons of matrix stuff: https://en.wikipedia.org/wiki/Rotation_matrix

    :param vec1: Vector one defined by [X Y Z]
    :param vec2: Vector two defined by [X Y Z]
    :return: 2D NumPy array of the rotation_matrix found
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    if s < 1e-5:
        # When vectors are almost parallel, align the cylinder with the direction vector
        rotation_matrix = np.eye(3)
    else:
        # Calculate the rotation matrix
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - np.dot(a, b)) / (s ** 2))
    return rotation_matrix


def find_smaller_z_value(graph: nx.Graph, node1: int, node2: int) -> np.array:
    """
    This function simply finds which of two nodes has a smaller Z value which is used in a fringe case of plotting.

    :param graph: NetworkX undirected graph which is stored as the design_graph object
    :param node1: 1st node to check
    :param node2: 2nd node to check
    """
    n1, n2 = graph.nodes[node1], graph.nodes[node2]
    P1z, P2z = n1['z'], n2['z']
    if P1z > P2z:
        # If P1 has a larger z value then we actually return n2 as the "bottom" point for the cylinder plot
        return np.array([n2['x'], n2['y'], n2['z']])
    else:
        return np.array([n1['x'], n1['y'], n1['z']])


def get_all_verts_in_graph(graph: nx.Graph, vertex_mapping: dict) -> np.array:
    """
    This function just returns all vertices found within a graph as a 2D NumPy array

    :param graph: NetworkX undirected graph which is stored as the design_graph object
    :param vertex_mapping: Dictionary of what vertex goes to which position in the graph.
    :return: NumPy array of vertices
    """
    all_verts = np.zeros((len(graph.nodes), 3))  # Create the array to store data
    for node in graph.nodes():
        pt = xyz_from_graph(graph=graph, node=node)
        all_verts[vertex_mapping[node]] = pt  # Breaking because this isn't technically correct logic. Node numbers always go up!
    return all_verts


def get_all_edges_in_graph(graph: nx.Graph) -> np.array:
    """
    This function just returns all edges found within a graph as a 2D NumPy array

    :param graph: NetworkX undirected graph which is stored as the design_graph object
    :return: NumPy array of all edges defined by the XYZ points [[[X1, Y1, Z1], [X2, Y2, Z2]], ...]
    """
    all_edges = []
    for edge in graph.edges():
        node1, node2 = edge[0], edge[1]
        pt1, pt2 = xyz_from_graph(graph=graph, node=node1), xyz_from_graph(graph=graph, node=node2)
        all_edges.append([pt1, pt2])
    return np.array(all_edges)


def unit_vector(P1: np.array, P2: np.array):
    """ Calculate the unit vector betewen two points """
    return (P1 - P2) / np.linalg.norm((P1 - P2))


def update_nodal_position(graph: nx.Graph, node: int, new_xyz: np.array, numDecimals: int) -> None:
    """
    This function will take a node and it's new XYZ positions and update the NetworkX graph to represent this.
    NOTE: This does NOT update the values in the triangular faces, another function controls this.

    :param graph: design_graph graph we are updating
    :param node: Node index #
    :param new_xyz: New XYZ Positions [0:X, 1:Y, 2:Z]
    :param numDecimals: Number of decimals to round values to
    :return: Nothing, just updates the design graph:
    """
    graph.nodes[node]['x'] = round(new_xyz[0], numDecimals)
    graph.nodes[node]['y'] = round(new_xyz[1], numDecimals)
    graph.nodes[node]['z'] = round(new_xyz[2], numDecimals)


def update_vertex_in_all_faces(all_faces: dict, node: int, new_xyz: np.array) -> None:
    """ This function updates the attribute values in stored MeshFace objects """
    for verts, face in all_faces.items():
        if node in list(verts):
            idx = verts.index(node)
            if idx == 0:
                access = 'v1'
            elif idx == 1:
                access = 'v2'
            else:
                access = 'v3'
            # If the node we are changing is in a face, we update the correct v1 / v2 / v3:
            setattr(face, access, new_xyz)



@jit(nopython=True)
def are_points_non_collinear(p1: np.array, p2: np.array, p3: np.array):
    """
    JIT compiled function that determines if 3 points are all collinear

    :param p1: [X Y Z] value of a point in 3D space
    :param p2: [X Y Z] value of a point in 3D space
    :param p3: [X Y Z] value of a point in 3D space
    :return: True if points are collinear, False if not
    """
    det = p1[0] * (p2[1] * p3[2] - p3[1] * p2[2]) - p2[0] * (p1[1] * p3[2] - p3[1] * p1[2]) + p3[0] * \
          (p1[1] * p2[2] - p2[1] * p1[2])
    return not np.isclose(det, 0)


def find_non_collinear_nodes(graph: nx.Graph, node: int, potential_nodes: list) -> tuple:
    """
    This function considers a graph, a node, and a list of potential nodes and determines which combination of the node
    and any two of the potential nodes results in a non-collinear set of three nodes. Otherwise, it return -1 -1 -1

    :param graph: design_graph graph we are updating
    :param node: Node index #
    :param potential_nodes: List of potential node indices to check to find a non-collinear set of three
    :return: Tuple of either 3 node indexes or (-1 -1 -1) signalling "no set found"
    """

    for combo in combinations(potential_nodes, 2):
        n1 = xyz_from_graph(graph=graph, node=node)
        n2 = xyz_from_graph(graph=graph, node=combo[0])
        n3 = xyz_from_graph(graph=graph, node=combo[1])
        if are_points_non_collinear(n1, n2, n3):
            return node, combo[0], combo[1]
    # Otherwise we can not create a plane:
    return -1, -1, -1


@jit(nopython=True)
def intersecting_edge_calc(P1: np.array, P2: np.array, P3: np.array, P4: np.array) -> np.array:
    """
    JIT compiled function that determines the distance between a set of points generated along the axes defined by
    P1 -> P2 and P3 -> P4. The goal of this function is to find the nearest distance between two lines in 3D space
    to determine if there is an intersection or not.

    :param P1: [X Y Z] value of a point in 3D space (used to defined axis P1 P2)
    :param P2: [X Y Z] value of a point in 3D space (used to defined axis P1 P2)
    :param P3: [X Y Z] value of a point in 3D space (used to defined axis P3 P4)
    :param P4: [X Y Z] value of a point in 3D space (used to defined axis P3 P4)
    :return: Numpy array of all distances calculated
    """
    axis1 = P2 - P1
    axis2 = P4 - P3

    # Calculate the step size for distributing points along the axis
    num_points = 15
    t_values = np.linspace(0, 1, num_points)
    # Calculate the points along axis1 and axis2:
    points_on_axis1 = P1 + t_values[:, np.newaxis] * axis1
    points_on_axis2 = P3 + t_values[:, np.newaxis] * axis2

    # Calculate distances from all points on axis1 to all points on axis 2 to create a 10x10 array for
    # collision box detection
    distances = np.sqrt(np.sum((points_on_axis1[:, np.newaxis, :] - points_on_axis2) ** 2, axis=2))
    return distances


def positive_value(x: float) -> bool:
    return x > 0  # Validate values across all dataclasses for all length (a, b, c) values


def sign(num):
    if num > 0:
        return 1
    elif num < 0:
        return -1
    else:
        return 0


#@jit(nopython=True)
def parallelepiped_volume(cell_type: str, a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> float:
    """
    This function simply calculates the volume of a parallelepiped but changes depending on the cell_type attribute as
    the volume can be more readily calculated through the equations below.
    """
    if cell_type == 'triclinic':
        magnitude = a * b * c
        angles = np.sqrt(1 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma) - (np.cos(alpha)**2) - (np.cos(beta)**2) - (np.cos(gamma)**2))
        volume = magnitude * angles
    elif cell_type == 'monoclinic':
        volume = a * b * c * np.sin(beta)
    elif cell_type == 'orthorhombic':
        volume = a*b*c
    elif cell_type == 'tetragonal':
        volume = a**2 * c
    elif cell_type == 'rhombohedral' or cell_type == 'trigonal':
        volume = a**3 * np.sin(alpha)
    elif cell_type == 'hexagonal':
        volume = a**2 * c * np.sin(alpha)
    elif cell_type == 'isometric' or cell_type == 'cubic':
        volume = a**3
    else:
        raise Exception('Invalid parallelepiped volume calculation')

    return volume


def find_graph_differences(graph1: nx.Graph, graph2: nx.Graph):
    """
    Function that determines any changes in the graph nodes / edges as well as any of the attributes between the nodes
    of the passed in graphs.

    :param graph1: design_graph before any design change was made
    :param graph2: design_graph after a grammar was applied to graph1
    :return: Set of differences in nodes, edges, and attributes
    """
    # Identify new nodes
    new_nodes = set(graph1.nodes) - set(graph2.nodes)

    # Identify new edges
    new_edges = set(graph1.edges) - set(graph2.edges)

    # Identify changes in node attributes
    node_attribute_changes = {}
    for node in graph1.nodes:
        if node in graph2.nodes:
            for attr in ['x', 'y', 'z']:
                if graph1.nodes[node].get(attr) != graph2.nodes[node].get(attr):
                    if node not in node_attribute_changes:
                        node_attribute_changes[node] = {}
                    node_attribute_changes[node][attr] = (graph1.nodes[node].get(attr), graph2.nodes[node].get(attr))

    return new_nodes, new_edges, node_attribute_changes
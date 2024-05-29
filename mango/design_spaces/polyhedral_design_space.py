"""
A.J. Vetturini
IDIG and MMBL
Carnegie Mellon University

This is the top level class responsible for maintaining the design features / qualities of a given design. This is
effectively an environment from which we are applying grammars (actions) to.
"""
# Import Modules
from dataclasses import dataclass
from typing import Optional, Union
from .bounding_box import TriclinicBox, MonoclinicBox, OrthorhombicBox, TetragonalBox, RhombohedralBox, HexagonalBox, CubicBox
from mango.mango_features.preserved_regions import PreservedEdge
from mango.mango_features.mesh_face import MeshFace
from mango.utils.mango_math import *
import trimesh
from itertools import combinations
from scipy.spatial import ConvexHull
from mango.utils.DNA_property_constants import BDNA
import pyvista as pv

@dataclass
class PolyhedralSpace(object):
    """
    This is a design space tasked with monitoring a polyhedral mesh for the export to automated design algorithms such
    as vHelix, DAEDALUS, and TALOS. This design space manages a 3D space via its design_graph

    Parameters
    ------------
    bounding_box : BoundingBox data type, created
    preserved : List of PreservedVertex or PreservedEdge objects
    excluded: An optional list of Sphere or RectangularPrism objects
    """
    bounding_box: Union[TriclinicBox, MonoclinicBox, OrthorhombicBox, TetragonalBox,
                        RhombohedralBox, HexagonalBox, CubicBox]
    preserved: list
    excluded: Optional[list]

    def __post_init__(self):
        self.design_graph = nx.empty_graph()  # Initialize this to an empty graph
        self.plotting_faces = {'preserved_region': [], 'other': []}
        self.nonterminal_nodes = []
        self.edge_divide_nodes = {}  # Tracking where all edge_divisions occur
        self.merge_label_dict = {}  # Used to track edge_divisions so we know when we can back out
        self.all_faces = {}  # Dictionary to tract MeshFaces thru stored tuples
        self.input_values = {}
        self.node_count = 0
        # This graph is used for the merge face grammar to guarantee that we can back out of certain topologies
        self.merge_label_graph = nx.empty_graph(create_using=nx.DiGraph)

        # Create the preserved regions by adding them to the design_graph:
        self.create_preserved_regions()

        # Tract any preserved edges such that we guarantee they are never removed:
        self.terminal_edges = []
        for i in self.preserved:
            if isinstance(i, PreservedEdge):
                p1 = node_from_xyz(graph=self.design_graph, xyz=i.v1)
                p2 = node_from_xyz(graph=self.design_graph, xyz=i.v2)
                self.terminal_edges.append((p1, p2))

        # Create a list of all initial vertices and edges used to validate the start shape:
        self.initial_nodes_pre_start = []
        for n in self.design_graph.nodes:
            self.initial_nodes_pre_start.append(n)

        self.mesh_file = None  # Initialize an object file for the mesh


    def increment_node_count(self) -> None:
        """
        Increments the node_count value by 1 so we can track the node values
        :return: Nothing, internally stored
        """
        self.node_count += 1


    def create_preserved_regions(self) -> None:
        """
        This function will use the passed in preserved region spaces and specify these regions in the design_graph

        :return: Nothing just updated the design_graph element
        """
        all_pts = {}
        for i in self.preserved:
            numNodes = i.preserved_region_as_graph.number_of_nodes()
            # We loop over the number of nodes in the region and change the relabel dict values as needed:
            for j in range(0, numNodes):
                pt = xyz_from_graph(graph=i.preserved_region_as_graph, node=j)
                pt = tuple(pt)
                if pt not in all_pts:
                    all_pts[pt] = self.node_count
                    i.relabel_dict[j] = self.node_count
                    self.increment_node_count()
                else:
                    # If the point already exists, we simply change the relabel dict to this value:
                    i.relabel_dict[j] = all_pts[pt]
            # Merging graphs
            new_preserved = nx.relabel_nodes(i.preserved_region_as_graph, i.relabel_dict)
            self.design_graph.add_nodes_from([(n, new_preserved.nodes[n]) for n in new_preserved.nodes])
            self.design_graph.add_edges_from(new_preserved.edges)


    def add_new_node(self, x: float, y: float, z: float, terminal: bool) -> None:
        """
        This function will create a new node in our design graph

        :param x: X position (nm)
        :param y: Y Position (nm)
        :param z: Z position (nm)
        :param terminal: Is node terminal (i.e., a node that can NOT be removed) or not
        :return: Nothing, modifies design_graph
        """
        self.design_graph.add_node(self.node_count, x=x, y=y, z=z, terminal=terminal)
        self.increment_node_count()


    def validate_start_shape(self, shape: Union[trimesh.Trimesh, pv.PolyData], preserved_edges: list,
                             og_verts: np.ndarray) -> tuple:
        """
        This function simply validates that the triangulated alphashape is able to be modified as expected by
        the application of various grammars

        :param shape: Mesh element generated by alphashape algorithm
        :param preserved_edges: A list of tuples that contain the nodes of any preserved edges.
        :param og_verts: Original list of vertices from the preserved_regions
        :return: (bool, list) : the bool returns if the shape is valid (true) or not (false) and the list contains
                                the list returns a list of tuples defining the graph nodes of the preserved regions
        """
        if not isinstance(shape, trimesh.Trimesh) and not isinstance(shape, pv.PolyData):
            raise Exception("An alphashape could not be generated for the provided input criterion and the optimization"
                            " process could not begin.")
        verts = shape.points
        tris = shape.faces.reshape(-1, 4)[:, 1:]
        # First we do a simple check to validate a mesh was triangulated:
        if len(tris) < 4:
            return True, -1  # If there are less than 4 triangles there is no way a 3D surface mesh is formed properly

        # First check if number of verices in mesh matches number of vertices passed in:
        if len(verts) != len(og_verts):
            return True, -1  # If we do not have the same number of points, then this mesh does not contain the input info

        # Now we must create a "map" of shape vertices. We do this because the delaunay triangulation might "move" the
        # indexed vertices around and so we want to validate we are checking the correct vertices:
        update_map = {}
        for idx, vert in enumerate(verts):
            curNode = node_from_xyz(graph=self.design_graph, xyz=vert)
            # If the current node in design_graph is not the same as the vert, we need to then update
            # start_shape_faces to conform to the proper value
            update_map[curNode] = idx

        # Now we create a loop that functionally does two things: 1) It validates that any preserved edges defined
        # by the user are in the starting mesh and 2) finds the new preserved_edge indices due to the "changed"
        # index positions from the map
        new_preserved_edges = []
        for i in preserved_edges:
            #n1, n2 = i[0], i[1]
            n1, n2 = update_map[i[0]], update_map[i[1]]
            new_preserved_edges.append((n1, n2))
            edge_found = False
            for f in tris:
                if n1 in f and n2 in f:
                    edge_found = True
                    break
            if not edge_found:
                return True, -1  # If we could not find a preserved edge we return True to re-arrange the mesh w new alpha

        # Otherwise, if we make it here, return False signalling "mesh passed checks, stop searching"
        return False, new_preserved_edges


    def triangulate_input_conditions(self) -> tuple:
        """
        This function uses the alphashape to create a start shape for the generative process based on the input
        conditions defined by the designer. Future work may want to consider different types of triangulations, but the
        alphashape generalizes well to a concave or convex set of input conditions and therefore I use it here.

        :return: A tuple of the faces of the triangulation and the related vertices of the mesh
        """
        # First gab the current design graph features:
        preserved, excluded = [], []
        for n in self.design_graph.nodes:
            # The current nodes will only contain preserved and excluded regions so just collect this info:
            pt = xyz_from_graph(graph=self.design_graph, node=n)
            if 'preserved_region' in self.design_graph.nodes[n]:
                preserved.append(pt)
            else:
                excluded.append(pt)  # Not currently used

        vert_list = np.array(preserved)

        # To generate a start shape we use the alpha-shape which requires an optimal alpha value which can be found:
        points = pv.wrap(vert_list)

        # Here we use an alpha shape to find a valid start shape within a while loop. Once a valid start shape is found
        # we break out. This is often just the convex hull of the entered points, but can expand. Future research
        # should look at this start shape condition as there are some assumptions made:

        keep_searching = True
        initial_alpha = 0.
        surface_mesh = None  # Init to get rid of annoying warning!
        while keep_searching:
            curSurf = points.delaunay_3d(alpha=initial_alpha)  # Note: reconstruct_surface() did not seem to work very well!
            surface_mesh = curSurf.extract_surface()
            # Validate the shape:
            keep_searching, new_edges = self.validate_start_shape(shape=surface_mesh, preserved_edges=self.terminal_edges,
                                                                  og_verts=vert_list)
            # If keep_searching is returned False then we update the terminal_edges:
            if not keep_searching:
                self.terminal_edges = new_edges
            initial_alpha += 0.1

            # Now we have a condition to check based on intiial_alpha:
            if initial_alpha >= 5:
                raise Exception('Unable to find a valid start condition for your input conditions, please reach out'
                                'to the developers to find a solution to your input case.')

        # After finding the correct surface_mesh:
        verts = surface_mesh.points
        tris = surface_mesh.faces.reshape(-1, 4)[:, 1:]

        return tris, verts


    def generate_start_shape(self) -> None:
        """
        This function adds all relevant edges on the initialized start_shape to the design_graph:
        """
        # The _ used to be a temp variable. I am currently leaving this as is just in case i need to debug later
        # since this isn't really computationally expensive. If i forget, OOPS!
        start_shape_faces, start_shape_vertices = self.triangulate_input_conditions()

        # Need to remap between start_shape vertices and the design_graph vertices so we create the correct
        # faces:
        update_map = {}
        for idx, vert in enumerate(start_shape_vertices):
            curNode = node_from_xyz(graph=self.design_graph, xyz=vert)
            # If the current node in design_graph is not the same as the vert, we need to then update
            # start_shape_faces to conform to the proper value
            update_map[curNode] = idx

        self.design_graph = nx.relabel_nodes(self.design_graph, update_map)  # Update positions properly per above map


        # First we extract all edges from the alpha shape and add them to the design_graph:
        for face in start_shape_faces:
            combs = list(combinations(face, 2))
            # First we add the edges to the design graph:
            for c in combs:
                edge = tuple(sorted(c))
                if edge not in self.design_graph.edges():  # IF THERE IS AN ISSUE CHECK HERE AJ, YOU GOT RID OF IF
                    self.design_graph.add_edge(edge[0], edge[1])
            # Now we add the faces to our plotting_faces list as well as the all_faces dictionary using the vertices:
            self.plotting_faces['other'].append(tuple(face))
            self.all_faces[tuple(face)] = MeshFace(
                v1=start_shape_vertices[face[0]],
                v2=start_shape_vertices[face[1]],
                v3=start_shape_vertices[face[2]], face_type='triangle'
            )

        # Now, we add the new faces to the merge_label_graph:
        # self.merge_label_graph.add_nodes_from(_)
        for node in start_shape_faces:  # Replaces above because apparently a TrackedArray worked but ndarray doesnt?
            self.merge_label_graph.add_node(tuple(node))


    def find_nonterminal_nodes(self) -> None:
        """
        This function will find all of the non-terminal nodes in the design and append to a list that will be used to
        apply a variety of production rules to manipulate geometry.

        Terminal nodes are those which belong to terminal geometries (preserved / excluded spaces). Nonterminal are
        all others

        :return: Nothing, just appends to an ever-changing list of nonterminal nodes
        """
        nonterminal_nodes = []
        for n in self.design_graph.nodes:
            # If a node is not marked as terminal and not already accounted for we go into this check
            if not self.design_graph.nodes[n]['terminal'] and n not in self.nonterminal_nodes:
                nonterminal_nodes.append(n)
        self.nonterminal_nodes = nonterminal_nodes  # Reset this list since nodes can be removed


    def create_trimesh_mesh(self) -> None:
        """
        This function is simply responsible for reconstructing a trimesh mesh file for the purpose of validating the
        mesh file. I deemed this "efficient enough" for the time being, but you could probably write some compiled
        code which does this much faster.
        """
        faces = np.array(list(self.all_faces.keys()), dtype=object)
        vertex_mapping = self.renumber_vertices(faces=faces)
        all_vertices = get_all_verts_in_graph(graph=self.design_graph, vertex_mapping=vertex_mapping)

        # Since nodes are constantly removed / added we need to re-number portions of the mesh:
        all_faces = []
        for row in faces:
            # We first use the mapping dictionary to ensure we are accessing the correct index positions:
            new_row = []
            for ind in list(row):
                new_row.append(vertex_mapping[ind])
            new_row_1 = [new_row[0], new_row[1], new_row[2]]
            all_faces.append(new_row_1)
        # Now with the faces re-mapped we can create the Trimesh mesh
        mesh = trimesh.Trimesh(vertices=all_vertices, faces=np.array(all_faces))
        mesh.fix_normals()

        # Assign it to the value of mesh:
        self.mesh_file = mesh


    def calculate_input_parameters(self, edge_lengths: list, routing_algorithm: str) -> dict:
        """
        This will return a standardized dictionary that will be used to allow for customized objective functions and
        design constraints, allowing for a more generalizable generative tool

        :param edge_lengths: A list of all edge lengths in units of nm
        :param routing_algorithm: Which algorithm is being used, options are currently: DAEDALUS and TALOS
        :return: dictionary {property : value} where property is the key value
        """
        # Calculate center-to-center distances of all edges:
        center_to_centers = []
        for e1, e2 in combinations(list(self.design_graph.edges()), 2):
            # We first need the centroid of each edge pair, so starting with e1:
            xyz1 = xyz_from_graph(graph=self.design_graph, node=e1[0])
            xyz2 = xyz_from_graph(graph=self.design_graph, node=e1[1])
            center_e1 = np.mean([xyz1, xyz2], axis=0)

            xyz3 = xyz_from_graph(graph=self.design_graph, node=e2[0])
            xyz4 = xyz_from_graph(graph=self.design_graph, node=e2[1])
            center_e2 = np.mean([xyz3, xyz4], axis=0)
            # Now, we want to find the distance between these centroids for d:
            dist = np.linalg.norm(center_e2 - center_e1)
            center_to_centers.append(dist)

        # Estimate scaffold length:
        total_length = np.sum(np.array(edge_lengths)) # Total edge length in nanometers
        if routing_algorithm == 'DAEDALUS':
            numBasepairs = 2 * int(total_length / BDNA.pitch_per_rise)  # We round to nearest integer
        elif routing_algorithm == 'TALOS':
            numBasepairs = 6 * int(total_length / BDNA.pitch_per_rise)  # We round to nearest integer
        else:
            raise Exception('Only valid routing algorithms supported are currently DAEDALUS and TALOS.')


        # Face Angles and Surface Areas:
        if self.mesh_file is None:
            self.create_trimesh_mesh()
        surface_areas = list(self.mesh_file.area_faces)

        ## Face Angles:
        face_angles = []
        for angles in self.mesh_file.face_angles:
            for i in angles:
                face_angles.append(np.rad2deg(i))

        # Principal axes based on the moment_inertia of the mesh file:
        inertia_axis_unit_vectors = self.mesh_file.principal_inertia_vectors

        # Convexity Measure:
        convex_hull = ConvexHull(self.mesh_file.vertices)

        # Calculate the volume of the convex hull
        convex_hull_volume = convex_hull.volume

        # Set values:
        self.input_values = {
            'edge_lengths': np.array(edge_lengths),
            'face_angles': np.array(face_angles),
            'estimated_scaffold_length': numBasepairs,
            'center_to_center': np.array(center_to_centers),
            'surface_area': np.array(surface_areas),
            'inertia_axis_unit_vectors': inertia_axis_unit_vectors,
            'volume': self.mesh_file.volume,
            'convex_hull_volume': convex_hull_volume,
            'design_graph': self.design_graph,
            'bounding_box': self.bounding_box
        }

        return self.input_values


    @staticmethod
    def renumber_vertices(faces: np.array) -> dict:
        """
        This function will take in the list of vertices and faces and re-number the faces to verify the index positions
        This is needed for writing out the .ply file to ensure proper file creation.

        :return: Dictionary containing the mapping of old to new vertices values
        """
        # First we need to get all the unique face indices:
        face_vertex_mapping = {}
        unique_faces = []
        for i in faces:
            for v in list(i):
                if v not in unique_faces:
                    unique_faces.append(v)
        # Next we sort this list and start "counting" through the list to create a mapping via a dictionary:
        unique_faces = sorted(unique_faces)
        vertex_counter = 0
        for i in unique_faces:
            face_vertex_mapping[i] = vertex_counter
            vertex_counter += 1
        return face_vertex_mapping


    def jsonify_design_space(self) -> dict:
        """
        This function will convert a design graph to a json-ified object such that I can plot it in the kodak
        visualization
        """
        # First call the trimesh mesh update for the current design state:
        self.create_trimesh_mesh()

        # Next loop over the trimesh mesh faces and vertices such that we can jsonify the object properly:
        vertices = self.mesh_file.vertices.tolist()
        faces = self.mesh_file.faces.tolist()

        # Create a dictionary containing vertices and faces
        mesh_data = {
            "vertices": vertices,
            "faces": faces
        }
        return mesh_data

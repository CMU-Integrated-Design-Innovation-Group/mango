"""
A.J. Vetturini
IDIG and MMBL
Carnegie Mellon University

This script contains the default (standard) design constraints used with the PolyhedralDesignSpace as well as the class
affording a designer the ability to create a custom defined constraint.
"""
from dataclasses import dataclass
from typing import Callable, Any, Dict, Tuple, Optional
import trimesh
from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
from mango.utils.mango_math import *
from mango.utils.DNA_property_constants import BDNA
from copy import deepcopy
import time
import pymeshfix as mf

@dataclass
class PolyhedralDefaultConstraints(object):
    """
    This default constraints class contains methods which can examine a design space and determine if the design is
    compliant w.r.t.: edge lengths, scaffold required to create the design, etc. If looking to constrain certain
    aspects of the design, please see the CustomDesignConstraint class.

    Parameters
    ------------
    extra_constraints : List of CustomDesignConstraint objects
    scaffold_routing_algorithm : Which automated algorithm is to be used on the exported PLY file, defualt DAEDALUS
    min_face_angle : Floating point minimal allowable face angle (in degrees). Too low of a value can cause issues
                     in the exported design.
    min_edge_length : Integer min number of nucleotides any one edge can have; default is 39 due to ATHENA restrictions.
    max_edge_length : Integer max number of nucleotides any one edge can have
    max_number_basepairs_in_scaffold : Integer value for max # of nucleotides to allow in a design. Note that here we
                                       estimate the # of nts in the scaffold, and there may be some error in exported
                                       designs.
    """
    extra_constraints: Optional[list] = None  # A list of CustomDesignConstraint objects that should also be checked
    scaffold_routing_algorithm: str = 'DAEDALUS'
    min_face_angle: float = 20  # Minimal face angle of a given polyhedral mesh in units of degrees
    min_edge_length: int = 39  # The default value presumes DAEDALUS whose min_edge_length is 39 basepair
    max_edge_length: int = 250  # User can change this if desired, I just picked a relatively random max value
    max_number_basepairs_in_scaffold: int = 7249  # Default value is M13

    def __post_init__(self):
        if self.scaffold_routing_algorithm not in ['DAEDALUS', 'TALOS']:
            raise Exception('Invalid scaffold routing algorithm passed in, the only supported options are currently '
                            'DAEDALUS and TALOS')

        if self.min_edge_length < 39:
            raise Exception('Invalid minimal edge length used with DAEDALUS scaffold routing algorithm, the minimal'
                            'value is 39 basepair.')

        # Create names:
        self.names = ['Outside Design Space', 'Vertex in Excluded', 'Edge in Excluded', 'Invalid Edge Length',
                      'Invalid Scaffold Length', 'Invalid Face Angle', 'Broken preserved edge', 'Intersecting edges',
                      'Intersecting faces']

        self.names_to_functions = {
            'Outside Design Space': (self.outside_design_space, ['design_space', 'all_edges']),
            'Vertex in Excluded': (self.vert_in_excluded_region, ['design_space', 'all_verts']),
            'Edge in Excluded': (self.edge_in_excluded_region, ['design_space']),
            'Invalid Edge Length': (self.invalid_edge_lengths, ['all_lengths']),
            'Invalid Scaffold Length': (self.used_too_much_scaffold, ['all_lengths']),
            'Invalid Face Angle': (self.invalid_face_angles, ['mesh']),
            'Broken preserved edge': (self.broken_preserved_edge, ['design_space']),
            'Intersecting edges': (self.find_intersecting_edges, ['design_space']),
            'Intersecting faces': (self.check_face_intersections, ['mesh'])
        }
        self.params = {}


    def update_params(self, design_space: PolyhedralSpace) -> None:
        """
        This function will update the proper value inside of self.params such that future design_constraints can be
        computationally efficiently checked

        :param design_space: The current design represented by its design_space to be validated for constraint
                             compliance
        """
        ## NOTE: This is a WIP, currently I am still literally re-calculating everything to get results rn:
        faces = np.array(list(design_space.all_faces.keys()), dtype=object)
        vertex_mapping = design_space.renumber_vertices(faces=faces)
        all_vertices = get_all_verts_in_graph(graph=design_space.design_graph, vertex_mapping=vertex_mapping)
        all_edge_lengths = calculate_design_edge_lengths(graph=design_space.design_graph)
        all_edges = get_all_edges_in_graph(graph=design_space.design_graph)
        design_space.create_trimesh_mesh()
        mesh = design_space.mesh_file
        self.params['design_space'] = design_space
        self.params['all_edges'] = all_edges
        self.params['all_verts'] = all_vertices
        self.params['all_lengths'] = all_edge_lengths
        self.params['mesh'] = mesh

    def create_constraint_dictionary(self, design_space: PolyhedralSpace) -> None:
        """
        This function creates the necessary dictionary which stores parameters (such as the trimesh) who are used
        in validating the design constraints of the problem

        :param design_space: The current design represented by its design_space to be validated for constraint
                             compliance
        """
        faces = np.array(list(design_space.all_faces.keys()), dtype=object)
        vertex_mapping = design_space.renumber_vertices(faces=faces)
        all_vertices = get_all_verts_in_graph(graph=design_space.design_graph, vertex_mapping=vertex_mapping)
        all_edge_lengths = calculate_design_edge_lengths(graph=design_space.design_graph)
        all_edges = get_all_edges_in_graph(graph=design_space.design_graph)
        design_space.create_trimesh_mesh()
        mesh = design_space.mesh_file
        self.params['design_space'] = design_space
        self.params['all_edges'] = all_edges
        self.params['all_verts'] = all_vertices
        self.params['all_lengths'] = all_edge_lengths
        self.params['mesh'] = mesh

    def check_constraints(self, design_space: PolyhedralSpace, constraint_set: list,
                          compare_space: PolyhedralSpace = None) -> tuple:
        """
        This is the top-level method which calls all constraints in the constraint_set (as well as any custom defined
        constraints).

        :param design_space: The current design represented by its design_space to be validated for constraint
                             compliance
        :param constraint_set: List of strings dictating which constraints are being used in the study
        :param compare_space: Design space of the prior design_space which is used to validate against the active design
        """
        # For all constraints we are checking we call the function with the params which is updated in the optimizer
        # code.
        store_time_dict = {}
        constraint_checks = []
        for constraint in constraint_set:
            function_to_call, param_names = self.names_to_functions[constraint]
            # Filter out only the parameters that the function requires
            valid_params = {param_name: self.params[param_name] for param_name in param_names}
            if constraint == 'Intersecting edges':
                valid_params['compare_space'] = compare_space  # Add the candidate design to speed up this calc
                constraint_failed, constraint_time = self.call_function_with_params(function_to_call, valid_params)
            else:
                constraint_failed, constraint_time = self.call_function_with_params(function_to_call, valid_params)
            store_time_dict[constraint] = constraint_time
            constraint_checks.append(constraint_failed)

        # Next we validate the custom defined constraints from the user:
        if self.extra_constraints:  # If we have constraints, then we will compute the parameters that a designer can
            # use for a custom constraint:
            input_params = design_space.calculate_input_parameters(edge_lengths=self.params['all_edges'],
                                                                   routing_algorithm=self.scaffold_routing_algorithm)
            for dc in self.extra_constraints:
                start = time.time()
                newCheck = dc.evaluate_constraint(input_parameters=input_params)
                constraint_time = time.time() - start
                constraint_checks.append(newCheck)
                store_time_dict[dc.name] = constraint_time
        # Finally, we return based on the values of the design constraints:
        if any(constraint_checks):
            # If any of these values are True we have an INVALID design and return False meaning the design is invalid:
            return False, store_time_dict
        else:
            # Otherwise we return True meaning GOOD (or non-design constraint breaking) design.
            return True, store_time_dict

    def call_function_with_params(self, function: Callable[..., Tuple[bool, float]], params: Dict[str, Any]) -> Tuple[
        bool, float]:
        return function(**params)

    # Individual Design Constraint Functions
    @staticmethod
    def outside_design_space(design_space: PolyhedralSpace, all_edges: np.array) -> tuple[bool, float]:
        """
        This method determines if a design is outside the bounding_box

        :param design_space: The current design represented by its design_space to be validated for constraint
                             compliance
        :param all_edges: List of all edge node pairs
        :return: True if the design is invalid, False if the design is valid. Also returns time of function execution
                 for internal logging / development purposes
        """
        start_time = time.time()
        box = design_space.bounding_box
        for edge in all_edges:
            # If we have an edge intersection then we return True signalling "invalid design"
            if box.shape.edge_intersecting_box(point1=edge[0], point2=edge[1]):
                constraint_time = time.time() - start_time
                return True, constraint_time

        # Otherwise we return False meaning the rule applied properly:
        constraint_time = time.time() - start_time
        return False, constraint_time

    @staticmethod
    def vert_in_excluded_region(design_space: PolyhedralSpace, all_verts: np.array) -> tuple[bool, float]:
        """
        This method determines if a vertex is inside of an excluded region which is not permissible.

        :param design_space: The current design represented by its design_space to be validated for constraint
                             compliance
        :param all_verts: List of all vertices in [X Y Z] as in [[X1 Y1 Z1], [X2 Y2 Z2], ... [XN YN ZN]]
        :return: True if the design is invalid, False if the design is valid. Also returns time of function execution
                 for internal logging / development purposes
        """
        start_time = time.time()
        # We must loop over all NP regions (and eventually the NP Regions)
        for r in design_space.excluded:
            # Now we loop over all the vertices to see if they are inside this region.
            for v in all_verts:
                check = r.point_inside_space(v)
                if check:
                    # If we are ever inside the NP Zone space, we just return True meaning INVALID DESIGN
                    constraint_time = time.time() - start_time
                    return True, constraint_time

        # If we get here, we return False meaning the rule applied properly:
        constraint_time = time.time() - start_time
        return False, constraint_time

    @staticmethod
    def edge_in_excluded_region(design_space: PolyhedralSpace) -> tuple[bool, float]:
        """
        This method determines if an edge crossed through a face in the excluded regions list.

        :param design_space: The current design represented by its design_space to be validated for constraint
                             compliance
        :return: True if the design is invalid, False if the design is valid. Also returns time of function execution
                 for internal logging / development purposes
        """
        # We must loop over all NP regions (and eventually the NP Regions)
        start_time = time.time()
        for r in design_space.excluded:
            # Now we loop over all the vertices to see if they are inside this region.
            for v1, v2 in design_space.design_graph.edges():
                # First we get the XYZ coordinates of each point:
                p1 = xyz_from_graph(graph=design_space.design_graph, node=v1)
                p2 = xyz_from_graph(graph=design_space.design_graph, node=v2)
                # Now we just pass to see:
                check = r.edge_intersecting_face(p1=p1, p2=p2)
                if check:
                    # If we are ever inside the NP Zone space, we just return True meaning INVALID DESIGN
                    constraint_time = time.time() - start_time
                    return True, constraint_time

        # If we get here, we return False meaning the rule applied properly:
        constraint_time = time.time() - start_time
        return False, constraint_time


    def invalid_edge_lengths(self, all_lengths: list) -> tuple[bool, float]:
        """
        This method simply verifies that all edge lengths in a design are between the min and max allowable value set
        by the user.

        :param all_lengths: List of edge length values [15, 20, 25] (nm)
        :return: True if the design is invalid, False if the design is valid. Also returns time of function execution
                 for internal logging / development purposes
        """
        start_time = time.time()
        # Reassign to a NumPy array for simple checking:
        all_lengths = np.array(all_lengths)
        # If any edge is below the min edge length threshold or above the max edge length threshold we return True
        ## First we need to convert from basepairs to nanometers for checking:
        if self.scaffold_routing_algorithm == 'DAEDALUS' or self.scaffold_routing_algorithm == 'TALOS':
            minEdge = BDNA.pitch_per_rise * self.min_edge_length
            maxEdge = BDNA.pitch_per_rise * self.max_edge_length
        else:
            minEdge = self.min_edge_length
            maxEdge = self.max_edge_length
        curMin, curMax = min(all_lengths), max(all_lengths)  # Min and max edge lengths in units of nm
        if curMin < minEdge or curMax > maxEdge:
            constraint_time = time.time() - start_time
            return True, constraint_time  # Return true meaning the design constraint failed
        # If we get here, we return False meaning the rule applied properly:
        constraint_time = time.time() - start_time
        return False, constraint_time


    def used_too_much_scaffold(self, all_lengths: list) -> tuple[bool, float]:
        """
        This method estimates the total number of nucleotides a scaffold would require to convert the design_space
        to a DNA origami design file using the specified automated routing algorithm (eg DAEDALUS). This method is
        simply an approximation and may over or under-shoot a design. Due to the stochastic nature of this algorithm,
        it would be prohibitively expensive to constantly compute the total number of nucleotides each design would
        require, whereas this is much more efficient.

        :param all_lengths: List of edge length values [15, 20, 25] (nm)
        :return: True if the design is invalid, False if the design is valid. Also returns time of function execution
                 for internal logging / development purposes
        """
        start_time = time.time()
        total_length = np.sum(np.array(all_lengths))
        total_length_in_number_of_nucleobases = int(total_length / BDNA.pitch_per_rise)

        # Depending on if we use DAEDALUS or TALOS we multiply this estimate by 2 or 6:
        if self.scaffold_routing_algorithm == 'DAEDALUS':
            total_length_in_number_of_nucleobases *= 2
        elif self.scaffold_routing_algorithm == 'TALOS':
            total_length_in_number_of_nucleobases *= 6

        if total_length_in_number_of_nucleobases > self.max_number_basepairs_in_scaffold:
            constraint_time = time.time() - start_time
            return True, constraint_time
        else:
            constraint_time = time.time() - start_time
            return False, constraint_time


    def invalid_face_angles(self, mesh: trimesh.Trimesh) -> tuple[bool, float]:
        """
        This method analyzes a mesh and assesses the faces of the mesh to assure none of the face angles are below the
        set minimal value.

        :param mesh: Trimesh mesh of the active design_space which auto-calculates all face angles
        :return: True if the design is invalid, False if the design is valid. Also returns time of function execution
                 for internal logging / development purposes
        """
        start_time = time.time()
        face_angles_degrees = np.rad2deg(mesh.face_angles)
        if np.min(face_angles_degrees) < self.min_face_angle:
            constraint_time = time.time() - start_time
            return True, constraint_time
        # If we get to the end, that means all face angles ARE valid so we return False:
        constraint_time = time.time() - start_time
        return False, constraint_time

    @staticmethod
    def broken_preserved_edge(design_space: PolyhedralSpace) -> tuple[bool, float]:
        """
        This method verifies that any PreservedEdge maintains as an edge during the generative process and is not
        manipulated via grammar application.

        :param design_space: The current design represented by its design_space to be validated for constraint
                             compliance
        :return: True if the design is invalid, False if the design is valid. Also returns time of function execution
                 for internal logging / development purposes
        """
        start_time = time.time()
        for terminal_edge in design_space.terminal_edges:
            dist, _ = length_and_direction_between_nodes(graph=design_space.design_graph, node1=terminal_edge[0],
                                                         node2=terminal_edge[1])
            shortest_paths = find_shortest_path(graph=design_space.design_graph, n1=terminal_edge[0],
                                                n2=terminal_edge[1])
            if min(shortest_paths) > dist:
                constraint_time = time.time() - start_time
                return True, constraint_time

        # If we make it here we signal False meaning the constraint passes (i know its the inverse of what you'd think)
        constraint_time = time.time() - start_time
        return False, constraint_time

    @staticmethod
    def find_intersecting_edges(design_space: PolyhedralSpace, compare_space: PolyhedralSpace = None) -> tuple[bool, float]:
        """
        This function instead considers the design as an invalid polyhedra where faces may intersect one another. This
        specific function will consider edges has cylindrical bounding boxes that may not collide with one another
        through a collision detection algorithm.

        :param design_space: The current design represented by its design_space to be validated for constraint
                             compliance
        :param compare_space: Design space of the prior design_space which is used to validate against the active design
        :return: True if the design is invalid, False if the design is valid. Also returns time of function execution
                 for internal logging / development purposes
        """
        start_time = time.time()
        nodes = list(design_space.design_graph.nodes())
        node_dict = all_node_xyz_dict(graph=design_space.design_graph, nodes=nodes)

        if compare_space is None:
            all_edges = list(design_space.design_graph.edges())
        else:
            # Here, we only calculate the distances of new edges. Furthermore, if a merge-rule occurred then we do
            # not need to bother with checking anything as it would be impossible for a collision to occur upon
            # a merge.
            new_nodes, new_edges, node_attribute_changes = find_graph_differences(design_space.design_graph,
                                                                                  compare_space.design_graph)

            all_edges = list(new_edges)  # Here we need to assign all_edges to be whatever the different edges are between
                              # compare_space and design_space

        edge_list_copy = deepcopy(list(design_space.design_graph.edges()))  # This will always be of all candidate
        # design edges

        # To start, we find all pairs of 2 edges that are not connected in the graph.
        for edge1 in all_edges:
            filtered_list = [tup for tup in edge_list_copy if edge1[0] not in tup]
            filtered_list = [tup for tup in filtered_list if edge1[1] not in tup]
            for edge2 in filtered_list:
                P1, P2 = node_dict[edge1[0]], node_dict[edge1[1]]
                P3, P4 = node_dict[edge2[0]], node_dict[edge2[1]]
                distances = intersecting_edge_calc(P1=P1, P2=P2, P3=P3, P4=P4)
                # Check if any distance is less than (2 * R (Diameter)) + e where e is some "interhelical" distance"
                # that is essentially a "safety" factor here of 0.25 nm.
                if np.any(distances < (BDNA.diameter + 0.25)):
                    # If at any point we find a collision, we return True meaning the design is invalid and stop
                    # the search
                    constraint_time = time.time() - start_time
                    return True, constraint_time  # Colliding
        # If all pairs are checked and not colliding, we return False meaning the design is good to go.
        constraint_time = time.time() - start_time
        return False, constraint_time


    @staticmethod
    def check_face_intersections(mesh: trimesh.Trimesh) -> tuple[bool, float]:
        """
        This function validates whether the faces of the mesh overlap or not signalling an "invalid" mesh.

        :param mesh: Trimesh mesh of the active design_space which auto-calculates all face angles
        :return: True if the design is invalid, False if the design is valid. Also returns time of function execution
                 for internal logging / development purposes
        """

        start_time = time.time()
        flag = False
        #mesh.show()
        tin = mf.PyTMesh()
        tin.load_array(mesh.vertices, mesh.faces)
        intersections = tin.select_intersecting_triangles()
        constraint_time = time.time() - start_time
        if len(intersections) > 0:
            flag = True

        return flag, constraint_time


@dataclass
class CustomDesignConstraint(object):
    """
    This class lets the user define a custom function which can be assigned as an extra design constraint. This is
    a slight WIP as I may want to convert this to a Python wrapper, but I am not sure about the logistics (or if this
    is even a good idea).

    The function can really be defined as anything, however it should be noted that there are no checks here to see how
    "good" a constraint is or is not, and that is at the discretion of the user.

    Parameters
    ------------
    name : String unique identifies of the constraint name
    design_constraint : A Python function defining the constraint. For example, if you have a constraint as:
                        def foo(input_vars, extra_params):
                            constraint_code

                        Then you would pass in CustomDesignConstraint(name='foo', design_constraint=foo,
                                                                      extra_params=extra_params)
    extra_params : A dictionary of constants / values that a user may want to use in the custom constraint
    """
    name: str  # Name of the custom design constraint from the user
    design_constraint: Callable = None  # This is the actual function defining the design constraint
    extra_params: dict = None  # Extra values that the user may pass in

    def evaluate_constraint(self, input_parameters):
        """ This function simply calls the Callable function """
        # If we are using a "simple" (or pre-programmed) function, we have to set it as a callable function:
        try:
            return self.design_constraint(input_parameters, self.extra_params)
        except TypeError:
            return self.design_constraint(input_parameters)

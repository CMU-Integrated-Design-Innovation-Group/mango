""" Test to validate that a design state can be represented in the caDNAno format in a precise study """
from mango.optimizers import MOSA, CustomDesignConstraint, PolyhedralDefaultConstraints, ObjectiveFunction
from mango.grammars import origami_grammars
from mango.features import Ramp, PreservedVertex
from mango.design_spaces import PolyhedralSpace, MonoclinicBox
import numpy as np
from itertools import combinations
from mango.utils.math import xyz_from_graph, length_and_direction_between_nodes
from scipy.spatial.distance import cdist
from copy import deepcopy



def variation_index(input_vars):
    edges_normalized = input_vars['edge_lengths'] / np.max(input_vars['edge_lengths'])
    if np.all(edges_normalized == edges_normalized[0]):
        standard_deviation = 0
    else:
        # Calculate the standard deviation
        #standard_deviation = np.std(edges_normalized) + np.std(angles_normalized)
        standard_deviation = np.std(edges_normalized)

    return standard_deviation


def cylinder_volume(graph, dna_diameter):
    total_volume = 0.
    for edge in graph.edges():
        cylinder_length, _ = length_and_direction_between_nodes(graph=graph, node1=edge[0], node2=edge[1])
        r_cyl_total = dna_diameter / 2  # Presume constant shell on all cylinders
        total_volume += (np.pi * r_cyl_total ** 2 * cylinder_length)
    return total_volume


def porosity_objective(input_vars, extra_params):
    cur_design = input_vars['design_graph']
    curBox = input_vars['bounding_box']
    total_volume = cylinder_volume(graph=cur_design, dna_diameter=extra_params['d'])
    curPorosity = curBox.shape.volume / total_volume
    return curPorosity


def cell_constraints(input_vars):
    bounding_box = input_vars['bounding_box']
    if bounding_box.a < 20 or bounding_box.a > 120:
        return True
    if bounding_box.b < 20 or bounding_box.b > 120:
        return True
    if bounding_box.c < 20 or bounding_box.c > 120:
        return True

    # Otherwise return False
    return False


def edge_cutoff_constraint(input_vars, extra_params):
    cur_design = input_vars['design_graph']
    threshold_to_check = extra_params['d'] + extra_params['threshold']

    def line_segment_points(start, end, threshold, num_points):
        """
        Generate points along the line segment from start + threshold to end - threshold.
        We do this to avoid comparing the small distance between two edges sharing a point.
        """
        direction_vector = end - start  # Vector from start to end
        length = np.linalg.norm(direction_vector)  # Length of the vector
        unit_vector = direction_vector / length  # Unit vector from start to end

        # Calculate the new starting and ending points, adjusted inwards by the threshold
        new_start = start + unit_vector * threshold
        new_end = end - unit_vector * threshold

        # Generate points between the adjusted start and end points
        return np.linspace(new_start, new_end, num_points)

    def disconnected_edge_pairs(nx_graph):
        for edge1, edge2 in combinations(nx_graph.edges(), 2):
            if not set(edge1).intersection(edge2):
                yield edge1, edge2


    for e1, e2 in disconnected_edge_pairs(cur_design):
        # Graph points e1 = [P1, P2] and e2 =[P3, P4]
        P1, P2 = xyz_from_graph(graph=cur_design, node=e1[0]), xyz_from_graph(graph=cur_design, node=e1[1])
        P3, P4 = xyz_from_graph(graph=cur_design, node=e2[0]), xyz_from_graph(graph=cur_design, node=e2[1])
        # Before comparing points, we also must consider if P1 is closer to P3 or P4 so that we are
        # fairly comparing distance arrays using cdist
        if np.linalg.norm(P3 - P1) > np.linalg.norm(P4 - P1):
            # If the distance to P4 from P1 is smaller than to P3, then we re-assign P3 and P4
            temp_value = deepcopy(P4)
            P4 = P3
            P3 = temp_value
        points1 = line_segment_points(P1, P2, threshold=5, num_points=5)
        points2 = line_segment_points(P3, P4, threshold=5, num_points=5)
        # Calculate all pairwise distances between points on the two line segments
        distances = cdist(points1, points2, 'euclidean')
        min_dist = np.min(distances)  # Find the smallest distance in the matrix
        # If any edge-to-edge distance is less than our cutoff distance, we reject the design:
        if min_dist < threshold_to_check:
            # If the minimal distance found is less than threshold, return True signalling "invalid design"
            return True

    # Otherwise, after checking all pairs, we return False signalling "valid design"
    return False


def create_box():
    X = 30
    Y = 40
    Z = 35
    random_seed = 48
    num_epochs = 75

    ## Define the bounding box, here we are targeting a specific value of "a" for a cubic box:
    new_box = MonoclinicBox(a=X, b=Y, c=Z, beta=80)

    # Define preserved regions at face midpoints:
    preserved_regions = []
    for midpoint in new_box.shape.midpoints:
        preserved_regions.append(PreservedVertex(v1=np.array(midpoint)))

    excluded_regions = []
    design_space = PolyhedralSpace(bounding_box=new_box, preserved=preserved_regions, excluded=excluded_regions,
                                   use_precise_geometry=False,
                                   #temp_geometry_path='./temp2',
                                   #automated_algorithm_software='/Users/kodak/Desktop/PERDIX-Mac/TALOS',
                                   #automated_name='TALOS'
                                   )

    # Next we import the grammar set and the default design constraints:
    grammar_set_1 = origami_grammars.TriangulationGrammars()
    grammar_set_2 = origami_grammars.ParallelepipedGrammars(cell_type=new_box.shape.cell_type)

    # Defining the custom constraints:
    # Lattice-constraints for the bounding box (min and max bound)
    cell_constraint = CustomDesignConstraint(name='cell constraints',
                                             design_constraint=cell_constraints)

    # Linear repulsive constraint to maintain edges a minimal distance apart from each other in an attempt to treat
    # the generated solutions as closely to the potential experimental configuration.
    r_cut_constraint = CustomDesignConstraint(name='Cutoff Distance Constraint',
                                              design_constraint=edge_cutoff_constraint,
                                              extra_params={'threshold': 3.0,
                                                            'd': 4.25})

    # Define out the constraint set:
    constraints = PolyhedralDefaultConstraints(min_face_angle=20,
                                               scaffold_routing_algorithm='TALOS',
                                               max_number_basepairs_in_scaffold=7249,
                                               extra_constraints=[cell_constraint,
                                                                  r_cut_constraint,
                                                                  ])

    # Create objective function objects:
    objective_1 = ObjectiveFunction(name='Porosity', objective_equation=porosity_objective, numDecimals=3,
                                    extra_params={'d': 6.})  # TALOS ~ 6nm diameter bundles
    objective_2 = ObjectiveFunction(name='Variation Index', objective_equation=variation_index, numDecimals=5)

    # Varied ramps:
    extension_ramp1 = Ramp(unique_id='Extension Ramp Triangulation', min_value=0.34, max_value=1.7,
                           max_number_epochs=num_epochs, min_steps_at_min=10)
    extension_ramp2 = Ramp(unique_id='Extension Ramp Parallelepiped', min_value=1, max_value=5,
                           max_number_epochs=num_epochs, min_steps_at_min=10)
    rotation_ramp = Ramp(unique_id='Rotation Ramp Both', min_value=1, max_value=10,
                         max_number_epochs=num_epochs, min_steps_at_min=10)
    rotation_ramp2 = Ramp(unique_id='Rotation Ramp Beta', min_value=1, max_value=5,
                          max_number_epochs=num_epochs, min_steps_at_min=10)

    # Create the MOSA optimizer with minimal input conditions to perform the test:

    ## DEFAULT CONSTRAINTS: (Not needed to be specified, just showing for clarity sakes)
    # Modified constraint set to accommodate the custom defined functions:
    cset = ['Outside Design Space', 'Vertex in Excluded', 'Edge in Excluded',
             'Invalid Edge Length', 'Invalid Scaffold Length', 'Invalid Face Angle',
             'Broken preserved edge', 'Intersecting edges', 'Intersecting faces']

    savename = 'TALOS_' + str(random_seed)
    optimizer = MOSA(
        design_space=design_space,
        grammars=[grammar_set_1, grammar_set_2],
        design_constraints=constraints,
        objective_functions=[objective_1, objective_2],
        SAVE_PATH="final_study",
        SAVE_NAME_NO_EXTENSION=savename,
        max_number_of_epochs=num_epochs,
        constraint_set=cset,
        extension_ramp={'Extend Vertex': extension_ramp1, 'Vary Box': extension_ramp2},
        rotation_ramp={'Edge Rotation': rotation_ramp, 'Rotate Box': rotation_ramp2},
        NT1=2000,
        NT2=2000,
        Na=1000,
        N_Bi=4000,
        random_seed=random_seed,
        max_time_of_optimization_minutes=600  # Just setting an upper bound in minute to terminate the process
    )

    # Start the process (which will automatically create the output file!)
    optimizer.begin_MOSA()


if __name__ == '__main__':
    create_box()
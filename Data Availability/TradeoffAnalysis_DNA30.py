"""
This script discusses how the multiobject studies were setup and run. This script can be run to reproduce the results
using the specified random seed #
"""
from mango.optimization_features.objective_function import ObjectiveFunction
from mango.design_spaces.bounding_box import MonoclinicBox
from mango.mango_features.preserved_regions import PreservedVertex
from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
from mango.grammars import origami_grammars
from mango.mango_features.grammar_ramp import Ramp
from mango.optimization_features import design_constraints
from mango.optimization_features.design_constraints import CustomDesignConstraint
from mango.optimizers.multiobjective_simulated_annealing import MOSA
import numpy as np
from itertools import combinations
from copy import deepcopy
from mango.utils.mango_math import xyz_from_graph
from scipy.spatial.distance import cdist

def objective1(input_vars):
    sa_to_v_ration = np.sum(input_vars['surface_area']) / input_vars['volume']
    return sa_to_v_ration

def objective2(input_vars):
    return input_vars['volume'] / input_vars['convex_hull_volume']

def cell_constraints(input_vars):
    bounding_box = input_vars['bounding_box']
    # If a is ever less than 20 we return True -> Failed constraint
    if bounding_box.a < 20 or bounding_box.a > 85:
        return True
    if bounding_box.b < 20 or bounding_box.b > 85:
        return True
    if bounding_box.c < 20 or bounding_box.c > 85:
        return True
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
    X = 50
    Y = 70
    Z = 55
    random_seed = 8
    num_epochs = 75

    ## Define the bounding box, here we are targeting a specific value of "a" for a cubic box:
    new_box = MonoclinicBox(a=X, b=Y, c=Z, beta=80)

    # Define preserved regions at face midpoints:
    preserved_regions = []
    for midpoint in new_box.shape.midpoints:
        preserved_regions.append(PreservedVertex(v1=np.array(midpoint)))

    excluded_regions = []
    design_space = PolyhedralSpace(bounding_box=new_box, preserved=preserved_regions, excluded=excluded_regions)

    # Next we import the grammar set and the default design constraints:
    grammar_set_1 = origami_grammars.TriangulationGrammars()
    grammar_set_2 = origami_grammars.ParallelepipedGrammars(cell_type=new_box.shape.cell_type)
    r_cut_constraint = design_constraints.CustomDesignConstraint(name='Cutoff Distance Constraint',
                                                                 design_constraint=edge_cutoff_constraint,
                                                                 extra_params={'threshold': 7,
                                                                               'd': 3.75})
    cell_constraint = CustomDesignConstraint(name='cell constraints', design_constraint=cell_constraints)

    constraints = design_constraints.PolyhedralDefaultConstraints(max_number_basepairs_in_scaffold=10000,
                                                                  min_face_angle=20,
                                                                  extra_constraints=[cell_constraint, r_cut_constraint])

    # Create objective function objects:
    objective_1 = ObjectiveFunction(name=f'Objective 1', objective_equation=objective1, numDecimals=3)
    objective_2 = ObjectiveFunction(name='Objective 2', objective_equation=objective2, numDecimals=3)

    # Varied ramps:
    extension_ramp1 = Ramp(unique_id='Extension Ramp Triangulation', min_value=0.34, max_value=3.4,
                           max_number_epochs=num_epochs, min_steps_at_min=10)
    extension_ramp2 = Ramp(unique_id='Extension Ramp Parallelepiped', min_value=0.4, max_value=2,
                           max_number_epochs=num_epochs, min_steps_at_min=10)
    rotation_ramp = Ramp(unique_id='Rotation Ramp Both', min_value=1, max_value=10,
                         max_number_epochs=num_epochs, min_steps_at_min=10)

    # Create the MOSA optimizer with minimal input conditions to perform the test:
    optimizer = MOSA(
        design_space=design_space,
        grammars=[grammar_set_1, grammar_set_2],
        design_constraints=constraints,
        objective_functions=[objective_1, objective_2],
        SAVE_PATH="./data_output/DNA30",
        SAVE_NAME_NO_EXTENSION="multiobjective_output2",
        max_number_of_epochs=num_epochs,
        extension_ramp={'Extend Vertex': extension_ramp1, 'Vary_a': extension_ramp2, 'Vary_b': extension_ramp2,
                        'Vary_c': extension_ramp2},
        rotation_ramp={'Edge Rotation': rotation_ramp, 'Rotate_beta': rotation_ramp},
        NT2=1000,
        Na=500,
        N_Bi=2000,
        random_seed=random_seed,
        max_time_of_optimization_minutes=300
    )

    # Start the process (which will automatically create the output file!)
    optimizer.begin_MOSA()


if __name__ == '__main__':
    create_box()





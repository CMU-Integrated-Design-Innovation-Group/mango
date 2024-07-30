""" Test to validate that a design state can be represented in the caDNAno format in a precise study """
from mango.optimizers.objective_function import ObjectiveFunction
from mango.design_spaces.bounding_box import MonoclinicBox
from mango.features.preserved_regions import PreservedVertex
from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
from mango.grammars import origami_grammars
from mango.features.grammar_ramp import Ramp
from mango.optimization_features import design_constraints
from mango.optimization_features.design_constraints import CustomDesignConstraint
from mango.optimizers.multiobjective_simulated_annealing import MOSA
import numpy as np
from mango.utils.math import length_and_direction_between_nodes
from itertools import combinations
from mango.utils.math import xyz_from_graph
from scipy.spatial.distance import cdist
from copy import deepcopy


def process_vhelix_rpoly_for_number_nucleotides(filepath):
    nucleotide_count = 0
    with open(filepath, 'r') as file:
        for line in file:
            # Strip the line of any leading/trailing whitespace
            stripped_line = line.strip()

            # If the line is empty and comments have ended, break the loop
            if stripped_line == "" and not comments:
                break

            # If the line is a comment (starts with '#'), skip it
            if stripped_line.startswith('#'):
                continue

            # If the line is not a comment and not empty, process it
            if stripped_line and not stripped_line.startswith('#'):
                comments = False
                parts = stripped_line.split()
                nucleotide_count += int(parts[2])
    # In vHelix's rpoly there are NO overhangs of any sort and every nucleotide is paired. The full design is twice the
    # nucleotide_count to account for the staples.
    return nucleotide_count

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

def molecular_weight_MSE(input_vars, extra_params):
    design_space = input_vars['PolyhedralSpace']
    if design_space.nucleotide_level_design is None:
        # If there is no cadnano design we have to manually assign a check to return a "very large" value. Because
        # we are minimizing in this study, returning a very large value will never be accepted as a valued state.
        return 1000000

    cadnano_design = design_space.nucleotide_level_design
    total_length_nts = 0
    for strand in cadnano_design.strands:
        for domain in strand.domains:
            total_length_nts += domain.dna_length()
    # For some computational savings, we will presume a constant 327 g/mol per nucleotide for DNA monophosphate
    # From https://www.thermofisher.com/us/en/home/references/ambion-tech-support/rna-tools-and-calculators/
    #                                   dna-and-rna-molecular-weights-and-conversions.html
    weight = total_length_nts * 327.0  # 327 in units of g/mol, this will give total MW of the structure:

    # Convert to MDa:
    weight /= 1e6

    # Calculate MSE to the target weight for optimization:
    error = np.abs((weight - extra_params['ideal_weight']) / extra_params['ideal_weight'])
    return error

def cell_constraints(input_vars):
    bounding_box = input_vars['bounding_box']
    # If a is ever less than 20 we return True -> Failed constraint
    if bounding_box.a < 20:
        return True
    if bounding_box.b < 20:
        return True
    if bounding_box.c < 20:
        return True
    return False

def scaffold_constraint(input_vars, extra_params):
    design_space = input_vars['PolyhedralSpace']
    if design_space.nucleotide_level_design is None:
        # If there is no cadnano design we return True automatically signalling "invalid design"
        return True

    cadnano_design = design_space.nucleotide_level_design
    total_length_nts = 0
    for domain in cadnano_design.scaffold.domains:
        total_length_nts += domain.dna_length()

    return total_length_nts > extra_params['max_scaffold_length']


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
    design_space = PolyhedralSpace(bounding_box=new_box, preserved=preserved_regions, excluded=excluded_regions,
                                   use_precise_geometry=True, temp_geometry_path='./test_temp/', automated_name='vHelix',
                                   # This must be specified upon execution
                                   automated_algorithm_software=r"C:\Users\Anthony\Desktop\vHelix_Test\bscor.bat")

    # Next we import the grammar set and the default design constraints:
    grammar_set_1 = origami_grammars.TriangulationGrammars()
    grammar_set_2 = origami_grammars.ParallelepipedGrammars(cell_type=new_box.shape.cell_type)

    # Defining the custom constraints:
    # Lattice-constraints for the bounding box (min and max bound)
    cell_constraint = CustomDesignConstraint(name='cell constraints',
                                             design_constraint=cell_constraints)
    # Maximal allowable scaffold length
    max_scaffold_constraint = CustomDesignConstraint(name='scaffold constraint',
                                                     design_constraint=scaffold_constraint,
                                                     extra_params={'max_scaffold_length':15000})

    # Linear repulsive constraint to maintain edges a minimal distance apart from each other in an attempt to treat
    # the generated solutions as closely to the potential experimental configuration.
    r_cut_constraint = design_constraints.CustomDesignConstraint(name='Cutoff Distance Constraint',
                                                                 design_constraint=edge_cutoff_constraint,
                                                                 extra_params={'threshold': 5.0,
                                                                               'd': 3.75})

    # Define out the constraint set:
    constraints = design_constraints.PolyhedralDefaultConstraints(min_face_angle=20,
                                                                  extra_constraints=[cell_constraint,
                                                                                     max_scaffold_constraint,
                                                                                     r_cut_constraint
                                                                                     ])

    # Create objective function objects:
    extra_params = {
        'd': 3.75,  # Diameter of helix bundles in design (DAED ~4, TALOS ~6)
    }
    objective_1 = ObjectiveFunction(name=f'Porosity Measure', objective_equation=porosity_objective, numDecimals=2,
                                    extra_params=extra_params)
    objective_2 = ObjectiveFunction(name='Molecular Weight', objective_equation=molecular_weight_MSE, numDecimals=4,
                                    extra_params={'ideal_weight': 15})  # Ideal weight in MDa units

    # Varied ramps:
    extension_ramp1 = Ramp(unique_id='Extension Ramp Triangulation', min_value=0.34, max_value=3.4,
                           max_number_epochs=num_epochs, min_steps_at_min=10)
    extension_ramp2 = Ramp(unique_id='Extension Ramp Parallelepiped', min_value=0.4, max_value=2,
                           max_number_epochs=num_epochs, min_steps_at_min=10)
    rotation_ramp = Ramp(unique_id='Rotation Ramp Both', min_value=1, max_value=5,
                         max_number_epochs=num_epochs, min_steps_at_min=10)

    # Create the MOSA optimizer with minimal input conditions to perform the test:
    ## NOTE: We are going to disable the invalid scaffold length and intersecting edges constraints which are on
    ##       by default. We are replacing these with custom constraints using the precise_geometry feature above.

    ## DEFAULT CONSTRAINTS:
    ## cset = ['Outside Design Space', 'Vertex in Excluded', 'Edge in Excluded',
    ##         'Invalid Edge Length', 'Invalid Scaffold Length', 'Invalid Face Angle',
    ##         'Broken preserved edge', 'Intersecting edges', 'Intersecting faces']

    # Modified constraint set to accommodate the custom defined functions:
    cset = ['Outside Design Space', 'Vertex in Excluded', 'Edge in Excluded',
            'Invalid Edge Length',  'Invalid Face Angle',
            'Broken preserved edge', 'Intersecting faces']
    optimizer = MOSA(
        design_space=design_space,
        grammars=[grammar_set_1, grammar_set_2],
        design_constraints=constraints,
        objective_functions=[objective_1, objective_2],
        SAVE_PATH="./data_output",
        SAVE_NAME_NO_EXTENSION="precise_geometry_example",
        max_number_of_epochs=num_epochs,
        constraint_set=cset,
        extension_ramp={'Extend Vertex': extension_ramp1, 'Vary_a': extension_ramp2, 'Vary_b': extension_ramp2,
                        'Vary_c': extension_ramp2},
        rotation_ramp={'Edge Rotation': rotation_ramp, 'Rotate_beta': rotation_ramp},
        NT2=2000,
        Na=1000,
        N_Bi=4000,
        random_seed=random_seed,
        max_time_of_optimization_minutes=360  # Just setting an upper bound in minute to terminate the process
    )

    # Start the process (which will automatically create the output file!)
    optimizer.begin_MOSA()


if __name__ == '__main__':
    create_box()
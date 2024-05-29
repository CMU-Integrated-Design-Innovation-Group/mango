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
from mango.utils.mango_math import length_and_direction_between_nodes

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

# Define objectives used in tradeoff:
def variation_index(input_vars):
    edges_normalized = input_vars['edge_lengths'] / np.max(input_vars['edge_lengths'])
    angles_normalized = input_vars['face_angles']
    if np.all(edges_normalized == edges_normalized[0]) and np.all(angles_normalized == angles_normalized[0]):
        standard_deviation = 0
    else:
        # Calculate the standard deviation
        # standard_deviation = np.std(edges_normalized) + np.std(angles_normalized)
        standard_deviation = np.std(edges_normalized)

    ## To prevent optimizer from getting stuck in a local region around a low index, we set the lower bound
    ## of "variation" to consider:
    return standard_deviation

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
    cell_constraint = CustomDesignConstraint(name='cell constraints', design_constraint=cell_constraints)
    constraints = design_constraints.PolyhedralDefaultConstraints(max_number_basepairs_in_scaffold=10000,
                                                                  min_face_angle=20,
                                                                  extra_constraints=[cell_constraint])

    # Create objective function objects:
    extra_params = {
        'd': 3.75,  # Diameter of helix bundles in design (DAED ~4, TALOS ~6)
    }
    objective_1 = ObjectiveFunction(name=f'Porosity Measure', objective_equation=porosity_objective, numDecimals=2,
                                    extra_params=extra_params)
    objective_2 = ObjectiveFunction(name='Variation Index', objective_equation=variation_index, numDecimals=4)

    # Varied ramps:
    extension_ramp1 = Ramp(unique_id='Extension Ramp Triangulation', min_value=0.34, max_value=3.4,
                           max_number_epochs=num_epochs, min_steps_at_min=10)
    extension_ramp2 = Ramp(unique_id='Extension Ramp Parallelepiped', min_value=0.4, max_value=2,
                           max_number_epochs=num_epochs, min_steps_at_min=10)
    rotation_ramp = Ramp(unique_id='Rotation Ramp Both', min_value=1, max_value=5,
                         max_number_epochs=num_epochs, min_steps_at_min=10)

    # Create the MOSA optimizer with minimal input conditions to perform the test:
    optimizer = MOSA(
        design_space=design_space,
        grammars=[grammar_set_1, grammar_set_2],
        design_constraints=constraints,
        objective_functions=[objective_1, objective_2],
        SAVE_PATH="./data_output",
        SAVE_NAME_NO_EXTENSION="multiobjective_output_Figure6",
        max_number_of_epochs=num_epochs,
        extension_ramp={'Extend Vertex': extension_ramp1, 'Vary_a': extension_ramp2, 'Vary_b': extension_ramp2,
                        'Vary_c': extension_ramp2},
        rotation_ramp={'Edge Rotation': rotation_ramp, 'Rotate_beta': rotation_ramp},
        NT2=2000,
        Na=1000,
        N_Bi=4000,
        random_seed=random_seed,
        max_time_of_optimization_minutes=300
    )

    # Start the process (which will automatically create the output file!)
    optimizer.begin_MOSA()


if __name__ == '__main__':
    create_box()





"""
This script uses multiple grammar sets as inputs to ShapeAnnealing to ensure it is able to discern between the sets
"""
# Import Modules
from mango.optimization_features.objective_function import ObjectiveFunction
from mango.design_spaces.bounding_box import TetragonalBox
from mango.mango_features.preserved_regions import PreservedVertex
from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
from mango.grammars.origami_grammars import TriangulationGrammars
from mango.optimization_features import design_constraints
from mango.mango_features.grammar_ramp import Ramp
import numpy as np
from mango.optimizers.single_objective_shape_annealing import ShapeAnneal


# Creating objective functions to use in the MOSA process:
def objective1(input_vars):
    sa_to_v_ration = np.sum(input_vars['surface_area']) / input_vars['volume']
    return sa_to_v_ration


if __name__ == '__main__':
    X = 50
    Y = 50
    Z = 70
    random_seed = 8
    num_epochs = 50  # Number of epochs to run study

    new_box = TetragonalBox(a=X, c=Z)
    preserved_regions = [
        PreservedVertex(v1=np.array([X / 2, 0, Z / 2])),
        PreservedVertex(v1=np.array([X / 2, Y, Z / 2])),
        PreservedVertex(v1=np.array([0, Y / 2, Z / 2])),
        PreservedVertex(v1=np.array([X, Y / 2, Z / 2])),
        PreservedVertex(v1=np.array([X / 2, Y / 2, 0])),
        PreservedVertex(v1=np.array([X / 2, Y / 2, Z])),
    ]

    excluded_regions = []

    design_space = PolyhedralSpace(bounding_box=new_box, preserved=preserved_regions, excluded=excluded_regions)

    # Next we import the grammar set and the default design constraints:
    grammar_set1 = TriangulationGrammars()
    constraints = design_constraints.PolyhedralDefaultConstraints()

    # Create objective function objects:
    objective_1 = ObjectiveFunction(name='SA : V Ratio', objective_equation=objective1)

    # Defining ramps:
    extension_ramp = Ramp(unique_id='Extension Ramp', min_value=0.34, max_value=6.8,
                          max_number_epochs=num_epochs, min_steps_at_min=10)

    rotation_ramp = Ramp(unique_id='Rotation Ramp', min_value=1, max_value=10,
                         max_number_epochs=num_epochs, min_steps_at_min=10)

    # Create the shape annealing optimizer with minimal input conditions to perform the test:
    optimizer = ShapeAnneal(
        design_space=design_space,
        grammars=grammar_set1,
        design_constraints=constraints,
        objective_function=objective_1,
        SAVE_PATH="./test_outputs",
        SAVE_NAME_NO_EXTENSION="ramp_test",
        extension_ramp={'Extend Vertex': extension_ramp},
        rotation_ramp={'Edge Rotation': rotation_ramp},
        random_walk_steps=100,
        max_number_of_epochs=num_epochs,
    )
    # Start the process (which will automatically create the output file!)
    optimizer.begin_annealing()


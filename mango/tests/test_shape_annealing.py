"""
This script performs a random walk of N steps to validate the grammars and design constraint capabilities. If we were
to get an error, there would be a flaw in the logic somewhere in the code if we consider testing across a large set of
N values.
"""
# Import Modules
from mango.optimization_features.objective_function import ObjectiveFunction
from mango.design_spaces.bounding_box import TetragonalBox
from mango.mango_features.preserved_regions import PreservedVertex
from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
from mango.grammars.origami_grammars import TriangulationGrammars
from mango.optimization_features import design_constraints
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
    grammar_set = TriangulationGrammars()
    constraints = design_constraints.PolyhedralDefaultConstraints()

    # Create objective function objects:
    objective_1 = ObjectiveFunction(name='SA : V Ratio', objective_equation=objective1)

    # Create the shape annealing optimizer with minimal input conditions to perform the test:
    optimizer = ShapeAnneal(
        design_space=design_space,
        grammars=grammar_set,
        design_constraints=constraints,
        objective_function=objective_1,
        SAVE_PATH="./test_outputs",
        SAVE_NAME_NO_EXTENSION="shape_annealing_test",
        random_walk_steps=100,
        max_number_of_epochs=10,
    )
    # Start the process (which will automatically create the output file!)
    optimizer.begin_annealing()


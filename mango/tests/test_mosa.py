"""
This script performs a random walk of N steps to validate the grammars and design constraint capabilities. If we were
to get an error, there would be a flaw in the logic somewhere in the code if we consider testing across a large set of
N values.
"""
# Import Modules
from mango.optimization_features.objective_function import ObjectiveFunction
from mango.design_spaces.bounding_box import TriclinicBox
from mango.mango_features.preserved_regions import PreservedVertex
from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
from mango.grammars import origami_grammars
from mango.optimization_features import design_constraints
from mango.optimizers.multiobjective_simulated_annealing import MOSA
import numpy as np


# Creating objective functions to use in the MOSA process:
def objective1(input_vars):
    sa_to_v_ration = np.sum(input_vars['surface_area']) / input_vars['volume']
    return sa_to_v_ration

def objective2(input_vars):
    return input_vars['volume'] / input_vars['convex_hull_volume']


if __name__ == '__main__':
    X = 50
    Y = 70
    Z = 90
    random_seed = 8
    num_epochs = 5

    new_box = TriclinicBox(a=X, b=Y, c=Z, alpha=70, beta=80, gamma=75)
    preserved_regions = []
    for midpoint in new_box.shape.midpoints:
        preserved_regions.append(PreservedVertex(v1=np.array(midpoint)))

    excluded_regions = []

    design_space = PolyhedralSpace(bounding_box=new_box, preserved=preserved_regions, excluded=excluded_regions)

    # Next we import the grammar set and the default design constraints:
    grammar_set = origami_grammars.TriangulationGrammars()
    constraints = design_constraints.PolyhedralDefaultConstraints(max_number_basepairs_in_scaffold=9999999)

    # Create objective function objects:
    objective_1 = ObjectiveFunction(name='SA : V Ratio', objective_equation=objective1, numDecimals=4)
    objective_2 = ObjectiveFunction(name='Convexity Ratio', objective_equation=objective2, numDecimals=4)

    # Create the MOSA optimizer with minimal input conditions to perform the test:
    optimizer = MOSA(
        design_space=design_space,
        grammars=grammar_set,
        design_constraints=constraints,
        objective_functions=[objective_1, objective_2],
        SAVE_PATH="./test_outputs",
        SAVE_NAME_NO_EXTENSION="MOSA_test",
        max_number_of_epochs=num_epochs,
        NT2=1000,
        random_seed=8,
        max_time_of_optimization_minutes=180
    )
    # Start the process (which will automatically create the output file!)
    optimizer.begin_MOSA()


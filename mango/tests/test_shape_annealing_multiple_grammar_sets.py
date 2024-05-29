"""
This script uses multiple grammar sets as inputs to ShapeAnnealing to ensure it is able to discern between the sets

Note:
    The hyperparameters for this test are set a bit "higher" resulting in a deeper search. If you are running this, be
    aware that this may take ~10-15 minutes depending on machine
"""
# Import Modules
from mango.optimization_features.objective_function import ObjectiveFunction
from mango.design_spaces.bounding_box import TetragonalBox
from mango.mango_features.preserved_regions import PreservedVertex
from mango.mango_features.excluded_regions import Sphere
from mango.mango_features.grammar_ramp import Ramp
from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
from mango.grammars.origami_grammars import TriangulationGrammars, ParallelepipedGrammars
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

    excluded_regions = [
        Sphere(diameter=15, center=np.array([X / 2, Y / 2, Z / 2]))
    ]

    design_space = PolyhedralSpace(bounding_box=new_box, preserved=preserved_regions, excluded=excluded_regions)

    # Next we import the grammar set and the default design constraints:
    grammar_set1 = TriangulationGrammars()
    grammar_set2 = ParallelepipedGrammars(cell_type='tetragonal')
    constraints = design_constraints.PolyhedralDefaultConstraints()

    # Create objective function objects:
    objective_1 = ObjectiveFunction(name='SA : V Ratio', objective_equation=objective1)

    # Varied ramps:
    extension_ramp1 = Ramp(unique_id='Extension Ramp Triangulation', min_value=0.34, max_value=6.8,
                           max_number_epochs=num_epochs, min_steps_at_min=10)
    extension_ramp2 = Ramp(unique_id='Extension Ramp Parallelepiped', min_value=0.4, max_value=3,
                           max_number_epochs=num_epochs, min_steps_at_min=10)

    rotation_ramp = Ramp(unique_id='Rotation Ramp Triangulation', min_value=1, max_value=5,
                         max_number_epochs=num_epochs, min_steps_at_min=10)



    # Create the shape annealing optimizer with minimal input conditions to perform the test:
    optimizer = ShapeAnneal(
        design_space=design_space,
        grammars=[grammar_set1, grammar_set2],
        design_constraints=constraints,
        objective_function=objective_1,
        SAVE_PATH="./test_outputs",
        SAVE_NAME_NO_EXTENSION="shape_annealing_multiple_grammar_sets",
        extension_ramp={'Extend Vertex': extension_ramp1, 'Vary_a': extension_ramp2, 'Vary_c': extension_ramp2},
        rotation_ramp={'Edge Rotation': rotation_ramp},
        random_walk_steps=1000,
        max_number_of_epochs=num_epochs,
        n=200,
        limit=75
    )
    # Start the process (which will automatically create the output file!)
    optimizer.begin_annealing()


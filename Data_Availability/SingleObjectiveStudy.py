"""
A.J. Vetturini
Sample use case of mango for generative design.
"""
from mango.design_spaces.bounding_box import CubicBox
from mango.mango_features.preserved_regions import PreservedEdge, PreservedVertex
from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
from mango.mango_features.grammar_ramp import Ramp
from mango.grammars.origami_grammars import TriangulationGrammars
from mango.optimization_features import design_constraints
from mango.optimizers.single_objective_shape_annealing import ShapeAnneal
import numpy as np
from mango.mango_features.excluded_regions import Sphere
from mango.optimization_features.objective_function import ObjectiveFunction
import multiprocessing
from mango.utils.mango_math import length_and_direction_between_nodes

def cylinder_volume(graph, dna_diameter):
    total_volume = 0.
    for edge in graph.edges():
        cylinder_length, _ = length_and_direction_between_nodes(graph=graph, node1=edge[0], node2=edge[1])
        r_cyl_total = dna_diameter / 2  # Presume constant shell "thickness" on all cylinders
        total_volume += (np.pi * r_cyl_total ** 2 * cylinder_length)
    return total_volume

def porosity_objective(input_vars, extra_params):
    cur_design = input_vars['design_graph']
    total_volume = cylinder_volume(graph=cur_design, dna_diameter=extra_params['d'])
    curPorosity = extra_params['cell_volume'] / total_volume
    return curPorosity


def edge_cutoff_constraint(input_vars, extra_params):
    cur_design = input_vars['center_to_center']
    # If any center-to-center edge is less than our cutoff distance (2.5 nm in this case), we reject the design:
    if min(cur_design) < extra_params['cutoff']:
        return True
    return False


def create_design_space(excluded_region_in_design: bool, savename: str, max_number_nucleotides: int, seed: int,
                        use_rcut: bool, cutoff_distance: float):
    random_seed = seed
    num_epochs = 100
    X = 50.

    new_box = CubicBox(a=X)
    preserved_regions = [
        PreservedVertex(v1=np.array([25, 0, 25])),
        PreservedVertex(v1=np.array([25, 50, 25])),
        PreservedVertex(v1=np.array([0, 25, 25])),
        PreservedVertex(v1=np.array([50, 25, 25])),
        PreservedEdge(v1=np.array([25., 35., 50]), v2=np.array([25., 15., 50])),
        PreservedEdge(v1=np.array([25., 35., 0.]), v2=np.array([25., 15., 0.])),
    ]

    if excluded_region_in_design:
        excluded_regions = [
            Sphere(diameter=10, center=np.array([25, 40, 25])),
            Sphere(diameter=10, center=np.array([25, 10, 25])),
        ]
    else:
        excluded_regions = []
    design_space = PolyhedralSpace(bounding_box=new_box, preserved=preserved_regions, excluded=excluded_regions)

    # Next we import the grammar set and the default design constraints:
    grammar_set = TriangulationGrammars()

    r_cut_constraint = design_constraints.CustomDesignConstraint(name='Cutoff Distance Constraint',
                                                                 extra_params={'cutoff': cutoff_distance},
                                                                 design_constraint=edge_cutoff_constraint)
    if use_rcut:
        constraints = design_constraints.PolyhedralDefaultConstraints(
            min_face_angle=18,
            max_number_basepairs_in_scaffold=max_number_nucleotides,
            extra_constraints=[r_cut_constraint]
        )
    else:
        constraints = design_constraints.PolyhedralDefaultConstraints(
            min_face_angle=18,
            max_number_basepairs_in_scaffold=max_number_nucleotides,
        )

    # Specify a ramp to search with:
    extension_ramp = Ramp(unique_id='Extension Ramp', min_value=0.34, max_value=3.4,
                          max_number_epochs=num_epochs, min_steps_at_min=10)

    rotation_ramp = Ramp(unique_id='Rotation Ramp', min_value=1, max_value=8,
                         max_number_epochs=num_epochs, min_steps_at_min=10)

    # Specify objective function:
    extra_params = {
        'd': 3.75,  # Diameter of helix bundles in design (DAED ~4)
        # Some presumed "additional" thickness added onto the diameter of d for a silicon shell if we presume coating of the design.
        'cell_volume': new_box.shape.volume  # This is held constant in this generative process
    }
    objective = ObjectiveFunction(name=f'Porosity Measure', objective_equation=porosity_objective,
                                  extra_params=extra_params)

    # Specify optimizer and return object:
    optimizer = ShapeAnneal(
        design_space=design_space,
        grammars=grammar_set,
        design_constraints=constraints,
        objective_function=objective,
        SAVE_PATH="./data_output",
        SAVE_NAME_NO_EXTENSION=savename,
        extension_ramp={'Extend Vertex': extension_ramp},
        rotation_ramp={'Edge Rotation': rotation_ramp},
        random_walk_steps=1000,
        max_number_of_epochs=num_epochs,
        n=1200,
        limit=500,
        random_seed=random_seed
    )
    return savename, optimizer



if __name__ == '__main__':

    # Create the output file (Used in Figure 4)
    # Re-run code above to regenerate the code we need i think?
    name, optimizer = create_design_space(excluded_region_in_design=True, savename='Single_Objective_Porosity',
                                          max_number_nucleotides=7249, seed=8, use_rcut=True, cutoff_distance=5.0)
    optimizer.begin_annealing()


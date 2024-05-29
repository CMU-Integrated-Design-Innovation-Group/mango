"""
A.J. Vetturini
Sample use case of mango for generative design.
"""
from mango.design_spaces.bounding_box import CubicBox
from mango.mango_features.preserved_regions import PreservedEdge, PreservedVertex
from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
from mango.grammars.origami_grammars import TriangulationGrammars
from mango.optimization_features import design_constraints
from mango.optimizers.single_objective_shape_annealing import ShapeAnneal
import numpy as np
from mango.mango_features.excluded_regions import Sphere
from mango.optimization_features.objective_function import ObjectiveFunction
import multiprocessing
from mango.utils.mango_math import length_and_direction_between_nodes
from itertools import combinations
from mango.utils.mango_math import xyz_from_graph
from scipy.spatial.distance import cdist
from copy import deepcopy

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


def create_design_space2(excluded_region_in_design: bool, savename: str, max_number_nucleotides: int, seed: int,
                         rcut_value: float):
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
    if rcut_value == -1:
        rcut_value = 0.01  # Set a minimal value that isn't truly zero just for simplicity sakes.
    r_cut_constraint = design_constraints.CustomDesignConstraint(name='Cutoff Distance Constraint',
                                                                 design_constraint=edge_cutoff_constraint,
                                                                 extra_params={'threshold': rcut_value,
                                                                               'd': 3.75})

    constraints = design_constraints.PolyhedralDefaultConstraints(
        min_face_angle=20,
        max_number_basepairs_in_scaffold=max_number_nucleotides,
        extra_constraints=[r_cut_constraint]
      )

    # Specify objective function:
    extra_params = {
        'd': 3.75,  # Diameter of helix bundles in design (DAED ~4, TALOS ~6)
        # Some presumed "additional" thickness added onto the diameter of d for a silicon shell if we presume coating of the design.
        'cell_volume': new_box.shape.volume  # This is held constant in this generative process
    }
    objective = ObjectiveFunction(name=f'Porosity Measure', objective_equation=porosity_objective,
                                  extra_params=extra_params)

    # Specify constraint set, optimizer and return object:
    cset = ['Outside Design Space', 'Vertex in Excluded', 'Edge in Excluded', 'Invalid Edge Length',
            'Invalid Scaffold Length', 'Invalid Face Angle', 'Broken preserved edge', #'Intersecting edges',
            'Intersecting faces'
            ]
    opt = ShapeAnneal(
        design_space=design_space,
        grammars=grammar_set,
        design_constraints=constraints,
        objective_function=objective,
        constraint_set=cset,
        SAVE_PATH="./data_output/InternalFeatures",
        SAVE_NAME_NO_EXTENSION=savename,
        extension_value_default=1.36,  # Make constant 4bp moves
        rotation_value_degrees_default=5,
        random_walk_steps=1000,
        max_number_of_epochs=num_epochs,
        n=400,
        limit=275,
        T_min=1e-24,  # Set this MUCH lower == deeper search but can be more ambiguous results wise
        max_time_of_optimization_minutes=420,
        random_seed=random_seed,
        use_convergence=False,
        print_progress=False,
    )
    return savename, opt



def start_process(annealer):
    annealer.begin_annealing()


if __name__ == '__main__':
    NUMBER_OF_CORES_TO_USE = 6  # This is 12 studies per 1 N_Studies
    max_scaffold = [5438, 7249, 10000]  # https://onlinelibrary.wiley.com/doi/10.1002/smll.201300701 for 5438 scaffold
    use_excluded_region = [True, False]
    use_rcut = ['None', '50']
    N_studies = 1
    seeds_for_studies = [18]
    run_these_studies = {}
    ct = 1

    for _ in range(N_studies):
        for scaffold_length in max_scaffold:
            for use_excluded in use_excluded_region:
                for rcut in use_rcut:
                    if rcut == '50':
                        temp = 5.0
                    else:
                        temp = -1
                    new_save_name = str(scaffold_length) + '_Excluded' + str(use_excluded) + '_Cutoff' + str(
                        rcut) + '_internal_study_' + str(ct)
                    name, optimizer = create_design_space2(excluded_region_in_design=use_excluded,
                                                           savename=new_save_name,
                                                           max_number_nucleotides=scaffold_length,
                                                           seed=seeds_for_studies[ct - 1],
                                                           rcut_value=temp)
                    run_these_studies[name] = optimizer
        ct += 1


    # Call the processing:
    annealer_objects = list(run_these_studies.values())

    with multiprocessing.Pool(processes=NUMBER_OF_CORES_TO_USE) as pool:
        pool.map(start_process, annealer_objects)

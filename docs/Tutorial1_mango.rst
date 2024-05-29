Tutorial 1: Setting up a generative design study using single objective optimization
====================================================================================

This tutorial is going to step through some of the basics on how the
mango generative design package is intended to be used and showcase what
it is capable of. The full series of three tutorials should give you
enough understanding of how the package is configured, taking you from
input of a bounding box through output of DNA origami design files.

This tutorial will focus on setting up the generative design space,
defining the optimization problem, configuring a custom constraint, and
finally running the generative design process.

Defining the inputs for the design space
----------------------------------------

The following blocks of code discuss how a design space is initialized
using only a minimal set of input conditions known as the bounding box,
a list of preserved regions, and an optional list of excluded regions.

The bounding box gives a bounding area to confine the space to in units
of nm x nm x nm The preserved regions are a list of PreservedVertex or
PreservedEdge objects that define material that is not to be removed
during the generative process. The excluded regions are optional, but
prevent DNA from being added to these areas (e.g., a nanoparticle)

.. code:: ipython3

    # Import the required packages:
    import numpy as np
    from mango.design_spaces.bounding_box import CubicBox
    from mango.mango_features.preserved_regions import PreservedEdge, PreservedVertex
    from mango.mango_features.excluded_regions import Sphere
    from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
    
    # Define a box that is a constant 50 x 50 x 50 nm:
    new_box = CubicBox(a=50)
        
    # Define preserved vertices at the mid-face of all 6 bounding box faces (units of nm)
    preserved_regions = [
            PreservedVertex(v1=np.array([25, 0, 25])),
            PreservedVertex(v1=np.array([25, 50, 25])),
            PreservedVertex(v1=np.array([0, 25, 25])),
            PreservedVertex(v1=np.array([50, 25, 25])),
            PreservedEdge(v1=np.array([25., 35., 50]), v2=np.array([25., 15., 50])),
            PreservedEdge(v1=np.array([25., 35., 0.]), v2=np.array([25., 15., 0.])),
        ]
    
    # We will prevent material from being added to a spherial region defined as below in units of nm.
    excluded_regions = [
                Sphere(diameter=10, center=np.array([25, 40, 25])),
                Sphere(diameter=10, center=np.array([25, 10, 25])),
            ]
    
    # Lastly, we pass in these values as the required input to the Polyhedral Design Space:
    design_space = PolyhedralSpace(bounding_box=new_box, preserved=preserved_regions, excluded=excluded_regions)

Defining a custom constraint
----------------------------

Before we setup the default constraints, we must first define a custom
constraint (if we have any). Here, I will us an example of a repulsive
force one might presume wireframe origami would exude in a given space.
Here, we are considering the closest point-to-point distance of all
non-connected bundles of DNA. If this distance is less than our
threshold, then we will reject that design as we may be presuming (in
this example) that these designs exceeding this threshold would not
physically form in solution. NOTE: This threshold is meant to serve as a
simplistic example and has no physical meaning

.. code:: ipython3

    # Import some default library modules:
    from itertools import combinations
    from scipy.spatial.distance import cdist
    from copy import deepcopy
    
    # Import mango utilities:
    from mango.utils.mango_math import xyz_from_graph, length_and_direction_between_nodes
    
    # We start by defining a function which takes in two values: input_vars and extra_params 
    # The value of input_vars is essentially re-calculated every design iteration, and the design_graph
    # is what we can use to determine if a design is valid or not. extra_params is effectively a list
    # of constants (but can also be a parameterized function that gets re-calculated)
    def edge_cutoff_constraint(input_vars, extra_params):
        cur_design = input_vars['design_graph']
        # We presume that the nearest distance is at least the diamater of BDNA + some threshold
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
    
        # For all two disconnected edges in the graph:
        ## This will presume that the minimal face angle constraint will keep the "connected" edges a
        ## reasonable distance apart.
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

.. code:: ipython3

    # Import design constraint features and assign:
    from mango.optimization_features import design_constraints
    
    # We define the custom constraint simple as:
    custom_constraint = design_constraints.CustomDesignConstraint(name='Cutoff Distance Constraint',
                                                                  design_constraint=edge_cutoff_constraint,
                                                                  extra_params={'threshold': 4.0, # units nm
                                                                                'd': 3.75})
    # NOTE: edge_cutoff_constraint must be written as is, do not use edge_cutoff_constraint()
    
    # Now we set up the default constraints and re-assign the min face and edge length:
    constraints = design_constraints.PolyhedralDefaultConstraints(
            min_face_angle=20,
            min_edge_length = 42,  # 42 basepairs for min length
            max_number_basepairs_in_scaffold=7249, 
        
            # Finally assign the custom constraint as a list of CustomDesignConstraint object(s):
            extra_constraints=[custom_constraint]
          )

Defining the objective of the optimization problem
--------------------------------------------------

Similarly to defining a custom design constraint, we must define out a
custom objective to be minimized. I am working on implementing the
maximizer code (if there is enough request), but also as a general note:
you can transform a maximization problem -> minimization problem by
taking the inverse of the maximizing function.

.. code:: ipython3

    from mango.optimization_features.objective_function import ObjectiveFunction
    # Function to estimate the volume of a cylinder-representing-DNA
    def cylinder_volume(graph, dna_diameter):
        total_volume = 0.
        for edge in graph.edges():
            cylinder_length, _ = length_and_direction_between_nodes(graph=graph, node1=edge[0], node2=edge[1])
            r_cyl_total = dna_diameter / 2  # Presume constant shell "thickness" on all cylinders
            total_volume += (np.pi * r_cyl_total ** 2 * cylinder_length)
        return total_volume
    
    # Our objective will simply divide the total bounding box volume by the volume of the DNA (and adding more DNA will lower this function!)
    def porosity_objective(input_vars, extra_params):
        cur_design = input_vars['design_graph']
        total_volume = cylinder_volume(graph=cur_design, dna_diameter=extra_params['d'])
        curPorosity = extra_params['cell_volume'] / total_volume
        return curPorosity
    
    # Specify objective function:
    extra_params = {
        'd': 3.75,  # Diameter of helix bundles in design (presume 2 helix bundle is about 4nm effective radius)
        'cell_volume': new_box.shape.volume  # This is held constant in this generative process
    }
    objective = ObjectiveFunction(name=f'Porosity Measure', objective_equation=porosity_objective,
                                  extra_params=extra_params)

Defining the optimizer
----------------------

Here I will discuss how the optimizer class is created. Note that the
only optimizers in this package (as of this writing) are simulated
annealing and multiobjective simulated annealing. Generally, the larger
the hyperparameters the “deeper” the search where the trade off is time
spent searching. However, we should note that this framework is truly
designed for conceptual / design exploration, and fine tuning (or
optimizing) a design with the currently developed grammars will likely
not be efficient.

I recommend starting simple and slowly “ramping up” the hyperparameters
to find a sweet spot of computation time and results analysis

.. code:: ipython3

    from mango.optimizers.single_objective_shape_annealing import ShapeAnneal
    from mango.grammars.origami_grammars import TriangulationGrammars
    
    opt = ShapeAnneal(
            design_space=design_space,
            grammars=TriangulationGrammars(),
            design_constraints=constraints,
            objective_function=objective,
            SAVE_PATH="./output_folder",
            SAVE_NAME_NO_EXTENSION='my_first_generated_design',
            extension_value_default=1.36,  # Make constant 4bp moves
            rotation_value_degrees_default=5,
            max_number_of_epochs=10, # This is very low, meant to run this in ~5/10 mins
            n=100,
            limit=25,  # Any more than 50% of moves leading to lower obj == COOL T!
            max_time_of_optimization_minutes=60,
            random_seed=8, # Stochastic algorithm means random seeds are important!
            print_progress=False,
        )

.. code:: ipython3

    # To start the generative process, you can simply run the following:
    opt.begin_annealing()


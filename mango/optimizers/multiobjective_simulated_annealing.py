"""
A.J. Vetturini
IDIG and MMBL
Carnegie Mellon University

Currently this optimizer is only available for the PolyhedralSpace; future work will potentially need to incorporate
different spaces and update this class to allow different spaces.

This script incorporates the MOSA algorithm as developed by Suppapitnam, Seffen, Parks, and Clarkson in "A Simulated
Annealing Algorithm for Multiobjective Optimization" (Engineering Optimization, 2000, Vol 33. pp. 59-85)

There are some modifications made to the algorithm in reference to Suppapitnarm, Parks, Shea, and Clarkson in
"Conceptual Design of Bicycle Frames by Multiobjective Shape Annealing".
"""

# Import Modules
from dataclasses import dataclass, field
from copy import deepcopy
import warnings
from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
from mango.grammars.origami_grammars import GrammarSet, CustomGrammarSet
from mango.optimization_features import design_constraints
from mango.utils.mango_math import *
import time
from typing import Union, List
import pandas as pd
from random import random, seed, randint
import os
import dill

@dataclass
class MOSA(object):
    """
    The optimizer is the "heart and soul" of the framework. It contains the logic for applying grammars and controlling
    the optimization-driven generative process. There are various hyperparameters that can be tuned which establish
    the depth of search and time spent searching.

    This data class contains the multiobjective shape annealing (MOSA) algorithm implementation and is set to solve a
    MOP with at least 2 objectives defined, but technically any N number of constraints can be input. However, MOSA has
    only been tested with 2 or 3 objectives, so any higher dimensionality must be analyzed with detail.


    Parameters
    ------------
    design_space : Initial PolyhedralSpace representing the design which is to be optimized
    grammars : A singular GrammarSet or a list of GrammarSets to use in the generative process
    design_constraints : The PolyhedralDesignConstraint object which constraint the optimizer
    objective_functions : List of at least 2 objective function objects to be minimized
    SAVE_PATH : Path of where to save the resultant generative .aj1 file
    SAVE_NAME_NO_EXTENSION : Filename of resultant .aj1 file
    constraint_set : List of constraint name strings that are checked by the optimizer by default

    Optional / Other input parameters that can be changed:
    ------------
    extension_value_default : Default node extension value to use (default: 0.34 nm)
    extension_ramp : Dictionary mapping a extension grammar to a Ramp element
    rotation_value_degrees_default : Default edge rotation value to use (default: 1 degrees)
    rotation_ramp: Dictionary mapping a rotation grammar to a Ramp element
    numDecimals : Number of decimal places to round to
    print_progress : Bool which determines if messages are printed to the terminal (default: TRUE)

    Multiobjective Simulated Annealing Hyperparameters that can be controlled
    Note that the default parameters here won't necessarily give a "good" result, they must be tuned!
    ------------
    NT1 : Number of random walk steps to initialize the temperature of the space (default: 2000)
    NT2 : Number of inner-loop mutations to consider before lowering temperature (default: 400)
    Na : Number of accepted design solutions to consider before lowering temperature (default: 400, 0.4*NT2)
    r_b : Return to base parameter in range (0, 1) and dictates frequency of changing active design (default: 0.95)
    acceptance_function : Controls how often "worse" moves are accepted (default: standard)
    cooling_schedule : Cooling schedule used during simulated annealing (default: huang / HRSV)
    N_Bi : Number of iterations before returning to base (default: 800, 2*NT2)
    N_Bi_LowerLimit : Lower limit on number of iterations before conducting return to base (default: 25)
    minimal_candidate_set_size : Lower limit on number of designs to consider when returning to base (default: 10)
    r_i : Parameter controlling how large the candidate set size is in range (0, 1) (default: 0.95)
    max_time_of_optimization_minutes : Max time spent generating a design (default: 60 minutes)
    T_min : Minimal temperature to end the annealing process (default: 1e-8)
    cooling_rate_geometric : Cooling rate to use if using geometric cooling_schedule (default: 0.9)
    delta_T : Triki annealing schedule adaptation rate, only used if using triki cooling_schedule (default: 0.8)
    """
    # Input parameters that must be passed in:
    design_space: PolyhedralSpace
    grammars: Union[GrammarSet, list]  # Pass a list of any grammar sets used.
    design_constraints: design_constraints
    objective_functions: list
    SAVE_PATH: str  # This is where the output dill file is saved to
    SAVE_NAME_NO_EXTENSION: str
    constraint_set: List[str] = field(
        default_factory=lambda: ['Outside Design Space', 'Vertex in Excluded', 'Edge in Excluded',
                                 'Invalid Edge Length', 'Invalid Scaffold Length', 'Invalid Face Angle',
                                 'Broken preserved edge', 'Intersecting edges', 'Intersecting faces'])

    # Various input parameters that can be changed:
    extension_value_default: float = 0.34  # Default value if ramp is not used
    extension_ramp: dict = field(default_factory=dict)
    rotation_value_degrees_default: float = 1
    rotation_ramp: dict = field(default_factory=dict)
    # Determines if a ramp is being used for rotation(True) which ignores above var.
    numDecimals: int = 5

    max_number_of_epochs: int = 50
    random_seed: int = None  # Random seed a user can specify

    ## Hyperparameters that can be tuned by user:
    NT1: int = 2000  # (1000 recommended for Initial temperature calculation as recommended by authors of MOSA)
    NT2: int = 500
    Na: int = 200  # Inner loop control, = 0.4NT2 per MOSA algo
    r_b: float = 0.95  # Return to base parameter in (0, 1) and dictates frequency of return
    acceptance_function: str = 'Standard'
    cooling_schedule: str = 'Huang'
    # (close to 1 == more returns over time == greater "depth")
    N_Bi: int = 1000  # 2 * NT2 is default. Controls initial # of returns. Larger value == deeper search
    N_Bi_LowerLimit: int = 25  # Lower limit for # of grammars to successfully apply before return to base
    minimal_candidate_set_size: int = 10
    r_i: float = 0.95  # Return fraction parameter, recommended at least 0.9 but can be in (0, 1)
    max_time_of_optimization_minutes: int = 60  # Number of minutes before exiting the optimization for protection
    print_progress: bool = True
    T_min: float = 1e-8  # Because of math the temperature will never truly be 0 so we need a stop point

    # Cooling Schedule Parameters depending on schedule used (Huang is recommended)
    cooling_rate_geometric: float = 0.8  # The cooling schedule for geometric (only used if Geoemtric schedule is used)
    delta_T: float = 0.8  # Used with the Triki annealing schedule and controls adaptation. May need to modify this.

    def __post_init__(self):
        ## Dummy variables that I initialize here to prevent users from "touching" them:
        self.total_checked = 0
        self.accepted_via_probability = 0
        self.sim_time_seconds = 0.0
        self.final_time = None
        self.obj_temp_dict = {}
        self.tracked_objectives_at_Ti = {}
        self.tracked_objectives_all_temp = {}
        self.temp_tracker = {}
        self.MOSA_archive = {}
        self.time_tracker_dict = {}
        self.archive_post_temperature_initialization = []
        self.pareto_animation = {}
        self.phi_r = 1.0  # Parameter that will be used across algorithm for return to base calcs, do not change from 1!
        self.objectives = []
        self.p_accept_list = []
        self.acceptance_tracker = []
        self.acceptance_per_epoch = []

        # Run initialization methods:
        if self.random_seed is None:
            self.random_seed = randint(1, 10000000)
            print(f'No random seed specified, using the following seed number: {self.random_seed}')
        # Initialize here:
        seed(self.random_seed)
        if isinstance(self.grammars, list):
            # In the case that multiple GrammarSets are prescribed, we need to "combine" them:
            # First validate they are GrammarSets:
            for g in self.grammars:
                if not isinstance(g, GrammarSet):
                    raise Exception(f'Each GrammarSet in the list of grammars must be a GrammarSet datatype, you '
                                    f'specified: {g}')
                # Also set the seed in this grammar_set
                g.set_seed(seed_number=self.random_seed)

            # If all are grammar sets then we use the "Combined" GrammarSet
            self.grammar_set = CustomGrammarSet(grammar_sets=self.grammars)
        elif isinstance(self.grammars, GrammarSet):
            # In the case a singular set of grammars are specifed (the most common case), we just set grammar_set to
            # grammars
            self.grammars.set_seed(seed_number=self.random_seed)
            self.grammar_set = self.grammars
        else:
            # If grammars is neither, then raise an error:
            raise ValueError(f"The value of grammars must be either a list of GrammarSet's or a singular GrammarSet, "
                             f"you specified: {self.grammars}")

        # Validate ramp conditions:
        if self.extension_ramp != {}:
            # If the user specifies an extension_ramp, then we must validate the values to ensure it is a proper
            # data structure:
            for grammar, ramp in self.extension_ramp.items():
                if grammar not in self.grammar_set.grammar_names:
                    raise Exception(
                        f'Invalid ramp specified for {grammar} as this grammar is not in the list of grammars '
                        f'specified to the problem.')
            # If we validate all ramps w/o error, then we just set the flag to true:
            self.use_extension_ramp = True
        else:
            self.use_extension_ramp = False

        # Repeat above for the rotation_ramp.
        if self.rotation_ramp != {}:
            for grammar, ramp in self.rotation_ramp.items():
                if grammar not in self.grammar_set.grammar_names:
                    raise Exception(
                        f'Invalid ramp specified for {grammar} as this grammar is not in the list of grammars '
                        f'specified to the problem.')
            self.use_rotation_ramp = True
        else:
            self.use_rotation_ramp = False
        self.initialize_objective_function_lists()
        self.verify_input_values()


    def update_design_space(self, new_design_space: PolyhedralSpace) -> None:
        """
        This function updates the active design_space object whenever return-to-base is enacted or whenever the space
        needs to be re-assigned
        """
        self.design_space = new_design_space


    def verify_input_values(self):
        """ Function that simply validates the input parameters before starting an optimization process. """
        if self.acceptance_function not in ['Standard', 'standard', 'logistic', 'Logistic', 'linear', 'Linear']:
            raise Exception('Invalid acceptance function input, only values are: standard, logistic, linear')

        if self.cooling_schedule not in ['geometric', 'Geometric', 'triki', 'Triki', 'Huang', 'huang', 'HRSV']:
            raise Exception('Invalid cooling schedule input, only values are: geometric, triki, huang')


    def initialize_objective_function_lists(self):
        """ Initializes lists and dictionaries for tracking relevent properties during optimization process """
        if len(self.objective_functions) < 2:
            raise Exception('Invalid list of objective functions used, the list must be at least length 2 (i.e., 2 '
                            'objectives must be defined for a multi-objective optimization problem!')

        for count, objective in enumerate(self.objective_functions):
            self.obj_temp_dict[count] = 100  # Initialize a temperature value of "100" as just a "dummy" value
            self.tracked_objectives_at_Ti[count] = []  # This is a list tracking the obj. func values for an epoch
            self.tracked_objectives_all_temp[count] = []
            self.temp_tracker[count] = []  # This list will track the temperatures at each epoch


    def initialize_time_tracker_dict(self):
        """ Initialize a dictionary to track time for internal development purposes """
        # first create keys for the grammars and constraints. We use "0" since a time-to-execute will never be less than
        # 0 and these are just tracking the execution-times of these functions.
        for grammar in self.grammar_set.grammar_names:
            self.time_tracker_dict[grammar] = 0

        for constraint in self.design_constraints.names:
            self.time_tracker_dict[constraint] = 0

        # Add in some defaults for the timer
        for f in self.objective_functions:
            nameMax = f.name + "_max_time"
            nameMin = f.name + "_min_time"
            # Initialize the max and min values to values that will get over-written during the tracking process
            self.time_tracker_dict[nameMax] = -1
            self.time_tracker_dict[nameMin] = 1e6


    def calculate_objective_function_values(self, design_space: PolyhedralSpace) -> list:
        """
        This function is responsible for evaluating the design_space using the input objective functions. It also
        records all execution times of the objective functions for internal development purposes.

        :param design_space: Design space that is currently being evaluated
        :returns: List of floating point values of the objective function valuations
        """
        all_edge_lengths = calculate_design_edge_lengths(graph=design_space.design_graph)
        input_params = design_space.calculate_input_parameters(edge_lengths=all_edge_lengths,
                                                               routing_algorithm=self.design_constraints.scaffold_routing_algorithm)
        obj_values = []
        for f in self.objective_functions:
            nameMax = f.name + "_max_time"
            nameMin = f.name + "_min_time"
            st_t = time.time()
            obj_val = f.evaluate_function(input_params)
            end_t = time.time()
            elapsed = end_t - st_t
            # RECORD TIMES
            if elapsed > self.time_tracker_dict[nameMax]:
                self.time_tracker_dict[nameMax] = elapsed
            if elapsed < self.time_tracker_dict[nameMin]:
                self.time_tracker_dict[nameMin] = elapsed
            # RECORD OBJECTIVE FUNCTION VALUE
            obj_values.append(obj_val)
        return obj_values


    @staticmethod
    def archive_datapoint(archive: list, test_point: tuple) -> tuple[bool, list]:
        """
        This function will take in the current archived points as well as a test point to determine if I need to archive
        the point (and do this archiving if needed).
        """

        # Internal function to determine if a point is dominating or not:
        def dominates(new_point, archived_point):
            # p1 dominates p2 if it's better in at least one objective and not worse in any objective
            # unpacking but in the future I should make this code more efficient, right now i have it for readability
            no_worse_in_any_objective = True  # Assume that the new_point is no worse in any objective to start
            strictly_better_in_at_least_one = False  # Assume False to start
            for f_new, f_old in zip(new_point, archived_point):
                if f_new < f_old:
                    strictly_better_in_at_least_one = True  # If a new point is better in at least one, we say True
                elif f_new > f_old:
                    no_worse_in_any_objective = False  # if f_new is ever larger than f_old in any objective, then
                    # it is no longer no worse in an objective

            # Now we just return these two booleans which determines dominance:
            return no_worse_in_any_objective and strictly_better_in_at_least_one

        ### CASE 1: If a candidate solution dominates any member(s) of the archive, those members are removed and the
        #           new solution is added.
        # Find dominated points in the archive and mark their indices for removal:
        to_remove = []
        for i, archive_point in enumerate(archive):
            if dominates(new_point=test_point, archived_point=archive_point):
                # This case means that the test point dominates the archive point meaning we need to remove archive
                # point from the pareto:
                to_remove.append(archive_point)
        # If to_remove has any points added then we know that our archive point dominates at least 1 point and is
        # therefore added to the archive (by returning True) and we remove the points in the to_remove list.
        if to_remove:
            return True, to_remove

        ### CASE 2: If a new solution is dominated by any members of the archive, it is not archived:
        for i, archive_point in enumerate(archive):
            if dominates(archive_point, test_point):
                # If any solution in the archive is dominating the test_point, we return False to not add to archive.
                return False, to_remove

        ### CASE 3: If a test_point does not dominate any solution in the archive, and it is not dominated by any
        #           member of the archive, we DO store the point, but we do NOT remove any from the list. However, here
        #           the to_remove will ALWAYS by of value [] (empty) so we can just return it like that:
        return True, to_remove


    def calculate_and_store_objective_function_values(self, design_space: PolyhedralSpace) -> None:
        """
        This function is called so that the objective functions are stored to analyze the optimization performance.

        :param design_space: Design space that is currently being evaluated
        """
        new_objectives = self.calculate_objective_function_values(design_space=design_space)
        for i in range(len(new_objectives)):
            self.tracked_objectives_at_Ti[i].append(new_objectives[i])
            self.tracked_objectives_all_temp[i].append(new_objectives[i])

        archive_data, remove_list = self.archive_datapoint(archive=list(self.MOSA_archive.keys()),
                                                           test_point=tuple(new_objectives))
        if archive_data:
            # This removal loop will only be needed if we are archiving a point, if we are not archiving a point we will
            # never be removing a point
            for remove_from_archive in remove_list:
                self.MOSA_archive.pop(remove_from_archive)

            # Archive is stored as (magnitudes of objective functions) : design_space object
            newKey = tuple(new_objectives)
            self.MOSA_archive[newKey] = deepcopy(design_space)

    def get_extend_and_rotation_values(self, epoch: int) -> tuple:
        """
        If a ramp is being used, this method gets the current value to use and return it

        :param epoch: Current epoch number to obtain extension value from ramp
        :return: Extension value and rotation value as a tuple: (extension, rotation)
        """
        extension_map = {}
        if self.use_extension_ramp:
            for grammar, ramp in self.extension_ramp.items():
                extension_value = ramp.current_ramp_value(epoch=epoch)  # We take in the max ramp value
                extension_map[grammar] = extension_value

        # Now, if the user is not using a ramp / for grammars not in the ramp we use the default value. Note that some
        # of these grammars dont even use an extension / rotation value but we return one anyways due to the structure
        # of the functions. The value is essentially "unused" but we still need to access it.
        for grammar in self.grammar_set.grammar_names:
            if grammar not in extension_map:
                extension_map[grammar] = self.extension_value_default

        # Repeat above for rotation_map
        rotation_map = {}
        if self.use_rotation_ramp:
            for grammar, ramp in self.rotation_ramp.items():
                extension_value = np.deg2rad(ramp.current_ramp_value(epoch=epoch))  # We take in the max ramp value
                rotation_map[grammar] = extension_value

        for grammar in self.grammar_set.grammar_names:
            if grammar not in rotation_map:
                rotation_map[grammar] = np.deg2rad(self.rotation_value_degrees_default)

        return extension_map, rotation_map


    def initialize_temperatures(self):
        """
        This function performs a random walk through the objective space to automatically assign the temperature used
        in the acceptance criterion for the first epoch.
        """
        # Start by evaluating the current start shape:
        self.calculate_and_store_objective_function_values(design_space=self.design_space)

        # In a random walk we will ALWAYS accept a move, but we need a counter and while loop since we can still
        # violate constraints
        accepted_moves = 0
        while_loop_counter = 0  # Prevent infinite loop

        # Get the extension value from the ramp. This will return the defaults if a ramp is not used as a note!
        extension_map, rotation_map = self.get_extend_and_rotation_values(epoch=0)

        # Random Walk:
        while accepted_moves < self.NT1:
            candidate_design_state = deepcopy(self.design_space)
            new_grammar = self.grammar_set.pick_random_grammar()
            if 'rotation' in new_grammar.lower() or 'rotate' in new_grammar.lower():
                ext = rotation_map[new_grammar]
            else:
                ext = extension_map[new_grammar]

            # Randomly apply grammar and validate the design constraints:
            grammar_applied_successfully, grammar_time = self.grammar_set.call_grammar_function(
                grammar_selected=new_grammar,
                design_space=candidate_design_state,
                extension_value=ext)
            # Update the design constraint values and calculate if design is compliant:
            if grammar_applied_successfully:
                # Update the design constraint values and calculate if design is compliant:
                self.design_constraints.update_params(design_space=candidate_design_state)
                passed_all_constraints, constraint_times = self.design_constraints.check_constraints(
                    design_space=candidate_design_state, constraint_set=self.constraint_set,
                    compare_space=self.design_space)
            else:
                # If the grammar didn't apply just set this to False. Doesn't need to be this way but ;-)
                passed_all_constraints = False

            if grammar_applied_successfully and passed_all_constraints:
                # If the grammar is applied successfully and the constraints are all passed then this new state replaces
                # the current design_space
                self.calculate_and_store_objective_function_values(design_space=candidate_design_state)
                self.update_design_space(new_design_space=candidate_design_state)
                accepted_moves += 1

            #  No matter what, we delete this candidate from memory so it re-instantiates in while loop for peace of
            #  mind regarding objects. If the design_space was accepted in if above, then this all works out!
            del candidate_design_state

            # For the sake of no infinite loops, I use a while loop counter here to ensure it doesn't get stuck:
            while_loop_counter += 1
            if while_loop_counter > (10 * self.NT1):  # Prevent infinite loop
                warnings.warn(f'Unable to apply {self.NT1} without violating design constraints which may lead to '
                              f'an inferior solution set, be aware!')
                accepted_moves = self.NT1 + 1  # Just break the loop condition...


        # After performing the random walk, we use the methodology from S. R. White from "Concepts of scale in
        #         # simulated annealing. IEEE International Conference of Computer Aided Design, Port Chester, New York.
        #         # pp. 646-651" to initialize each temperature of each objective function:
        count = 0
        for key, val in self.tracked_objectives_at_Ti.items():
            newTemp = round(float(np.std(val)), self.numDecimals)
            self.obj_temp_dict[key] = newTemp
            self.tracked_objectives_at_Ti[key] = []  # Reset these lists as we will return to base after initializing
            # Check and make sure the temperature is not 0:
            if np.isclose(self.obj_temp_dict[key], 0):
                obj = self.objective_functions[count]
                raise Exception(f'MOSA could not find any variation in the objective function valuations for the'
                                f'specified function {obj.name}. Please increase NT1 or verify the objective function '
                                f'is actually evaluating at different values.')
            count += 1

        # We store the post-random walk Pareto for visualization purposes:
        self.archive_post_temperature_initialization = list(self.MOSA_archive.keys())


    def calculate_pareto_isolation(self) -> pd.DataFrame:
        """
        This function is used to rank and calculate the isolations of the pareto front during optimization.
        :return: Pandas dataframe containing candidate designs to return to base to
        """
        pareto = list(self.MOSA_archive.keys())
        # First I am going to split the pareto into lists for each objective function:
        transposed_tuples = list(zip(*pareto))
        maxBRF = [max(column) for column in transposed_tuples]
        minBRF = [min(column) for column in transposed_tuples]


        # Loop over all solutions
        isolation_dict = {}  # Dictionary to track measures of isolation
        for i in pareto:
            # Our inner loop is for all other solutions in the pareto
            Iij = 0
            for j in pareto:
                if i == j:
                    # Here nothing happens because we don't compare the isolation of the same point
                    continue
                else:
                    # Calculate isolation value for this index:
                    for f1, f2, FMAX, FMIN in zip(i, j, maxBRF, minBRF):
                        if FMAX == FMIN:
                            pass
                            # If FMAX and FMIN are the same then there is only one point in the Pareto set so we prevent
                            # a divide-by-zero error by just passing along.
                        else:
                            Iij += (((f1 - f2) / (FMAX - FMIN))**2)

            # Before continuing, we store this isolation measure:
            isolation_dict[i] = Iij

        # After calculating all isolation values, we then rank in order of decreasing isolation distance WITH THE
        # EXCEPTION of the extreme solutions (i.e. solutions corresponding to extrema in the trade-off / solutions
        # with the lowest objective function. Sorting the dictionary:
        points, isolations = [], []
        for k, v in isolation_dict.items():
            # Since I am currently using a relatively simple objective function space, I will not remove the extrema
            # datapoints, as these do not explicitly constitute a barely feasible solution.
            if any(val in minBRF for val in k):
                # If the point we are checking contains a point in the extrema (which in this case is a MINIMIZATION
                # problem), then we are not going to choose it in our return to base
                continue
            else:
                # If the point is not an extrema, we are going to add it to this 2D list:
                points.append(k)
                isolations.append(v)
        # For safety and if NT1 is too small, we check to see if points and isolations are empty, and if so, we add
        # the isolation dict values:
        if points == [] and isolations == []:
            for k, v in isolation_dict.items():
                points.append(k)
                isolations.append(v)

        # Combining into a dataframe to sort by descending order per MOSA algorithm:
        df = pd.DataFrame({'points': points, 'isolations': isolations})
        df = df.sort_values(by='isolations', ascending=False)
        datapoints = df.shape[0]
        selection_set_size = round(self.phi_r * datapoints)  # Rounding because we need an integer # of points to pick
        # We also incorporate a minimal selection set size to prevent ever focusing in on just a few solutions:
        if selection_set_size < self.minimal_candidate_set_size:
            selection_set_size = self.minimal_candidate_set_size

        # We select the first "selection_set_size" number of points to be into "candidate solutions" we will choose from
        candidate_list = df[:selection_set_size]

        # After creating the candidate list, we update values for future selections. But we only do this when we are not
        # performing the initial walk through space
        self.phi_r *= self.r_i

        return candidate_list


    def return_to_base(self) -> None:
        """
        The return to base functionality is used to expose the trade-off between objective functions. The return to base
        will randomly change the active design state to somewhere else within the Pareto so that we optimize the
        pareto front

        Here we utilize the "intelligent" return-to-base strategy laid out in Conceptual Design of Bicycle Frames by
        Multiobjective Shape Annealing. The purpose is to prefer to "return to base" to extreme solutions (or highly
        isolated solutions) to explore around these unknown areas of the design space
        """
        # First we calculate the pareto isolation parameter to determine where we will "return to base" to search
        candidate_list = self.calculate_pareto_isolation()

        # Finally, we randomly select a single candidate to become the new active state:
        new_state = candidate_list.sample(n=1)
        new_state = tuple(new_state['points'])  # Convert the datatype...
        store_state = tuple(new_state[0])
        new_state_data_object = self.MOSA_archive[store_state]

        # Lastly, before returning the data object, we update the N_Bi term which dictates how often we will return to
        # base during the search. In general, the longer the search the more often we return to base to exploit the
        # tradeoffs of the objective functions:

        # We only lower these values once we start the actual MOSA search
        self.N_Bi = int(round(self.r_b * self.N_Bi))
        if self.N_Bi < self.N_Bi_LowerLimit:
            self.N_Bi = self.N_Bi_LowerLimit  # If we ever go below the lower limit, we just set N_Bi to the limit

        self.update_design_space(new_design_space=new_state_data_object)


    def keep_annealing(self, start_time: float, numEpochs: int) -> bool:
        """
        This function determines when the MOSA algorithm should end based on either a set number of epochs, a length
        of time, or if the temperature hyperparameters get too low

        :param numEpochs: Current epoch of MOSA
        :param start_time: Start time of how long the inner while loop has been running
        :return: True = Keep running, False = END
        """
        elapsed_time = time.time() - start_time
        hours, remainder = divmod(elapsed_time, 3600)
        minutes, seconds = divmod(remainder, 60)
        time_formatted = "{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds)

        ### First check if temperatures are low enough to end the run:
        for temp in self.obj_temp_dict.values():
            if temp < self.T_min:
                if self.print_progress:
                    print(f'MOSA ended due to temperatures dropping below threshold of {self.T_min}')
                self.final_time = time_formatted
                return False

        ### Next check if max time has elapsed
        if (time.time() - start_time) >= (self.max_time_of_optimization_minutes * 60):
            # If the current elapsed time is greater than the threshold in seconds (why i multiply by 60 above), we end
            if self.print_progress:
                print(f'MOSA ended due to max time input of {self.max_time_of_optimization_minutes} minutes')
            self.final_time = time_formatted
            return False

        if numEpochs > self.max_number_of_epochs:
            if self.print_progress:
                print(f'MOSA ended due to max number of epochs of {self.max_number_of_epochs} being met')
            self.final_time = time_formatted
            return False
        else:
            # Otherwise we return True meaning we keep looping
            return True


    def update_constraint_time_tracker(self, new_times: dict) -> None:
        """ Updates time tracked with new_times dictionary """
        for dc_name, ti in new_times.items():
            try:
                if ti > self.time_tracker_dict[dc_name]:
                    self.time_tracker_dict[dc_name] = ti
            except KeyError:
                self.time_tracker_dict[dc_name] = ti


    def attempt_to_archive(self, candidate_design: PolyhedralSpace) -> bool:
        """
        This function will be used to determine if a point is archived or not during the inner while loop of MOSA

        :param candidate_design: Design space after applying a grammar that is being considered.
        :return: True if the point was ACCEPTED and False if the point was REJECTED
        """
        # Calculate objective function values and pass to archive function:
        objective_vals = self.calculate_objective_function_values(design_space=candidate_design)
        archive_data, remove_list = self.archive_datapoint(archive=list(self.MOSA_archive.keys()),
                                                           test_point=tuple(objective_vals))

        # If we are archiving the data point we use the same logic to update the dictionary depending on dominance
        if archive_data:
            # First store the objective functions values into the proper tracker lists:
            for i in range(len(objective_vals)):
                self.tracked_objectives_at_Ti[i].append(objective_vals[i])
                self.tracked_objectives_all_temp[i].append(objective_vals[i])

            # This removal loop will only be needed if we are archiving a point, if we are not archiving a point we will
            # never be removing a point as remove_list will be blank.
            for remove_from_archive in remove_list:
                self.MOSA_archive.pop(remove_from_archive)

            # Archive is stored as (magnitudes of objective functions) : design_space object
            self.MOSA_archive[tuple(objective_vals)] = deepcopy(candidate_design)

            # Finally, we return True as if a datapoint is archived, then that state is always accepted per MOSA:
            return True

        # If we do not archive, we return False to do a probability test:
        return False


    def p_accept(self, new: float, old: float, temperature: float) -> float:
        """
        This function will calculate the individual probability for MOSA

        :param new: Objective function value for objective i for the MODIFIED design state
        :param old: Objective function value for objective i for the ORIGINAL design state
        :param temperature: Current temperature for objective i
        """
        # Use this to catch computational errors in the exponential:
        overflow_check = ((new - old) / temperature)  # Just checking this value to make sure it's calculating out
        if self.acceptance_function == 'standard' or self.acceptance_function == 'Standard':
            # I first check to see if new - old / temperature is greater than a threshold as we are passing these to a
            # exponential function. This is because I want to avoid the numpy "very large value" warnings as well as
            # because i want to plot these probability values. A very large value for temp means that the move we are
            # evaluating is VERY positive to the objective function we are testing (it's temperature).

            if overflow_check <= 4.6:
                # This threshold of 4.6 correspond to about ~100 (e^-(new-old/temperature) = e^-4. = 100)
                # I will return "100" as that means we are going to accept this state as the new active state
                return 1
            else:
                return (-1 * (new - old)) / temperature

        elif self.acceptance_function == 'logistic' or self.acceptance_function == 'Logistic':
            if overflow_check <= 4.6:
                # This threshold of 4.6 correspond to very large values of e, and the larger this e value is the actual
                # value of temp used in Logistic APF goes to 0 since it uses e(x) instead of e(-x)
                temp = 0
            else:
                temp = np.exp((new - old) / temperature)  # Logistic Curve Acceptance Probabilty Function
            return 2 / (1 + temp)

        elif self.acceptance_function == 'linear' or self.acceptance_function == 'Linear':
            if overflow_check <= 4.6:
                return 1  # If the threshold is met, we return the value 1 since the APF is min([1, exp])
            else:
                temp = np.exp(((-1 * (new - old)) / temperature))  # Scalar linear Acceptance Probabilty Function
                return min([1, temp])
        else:
            raise Exception('Invalid acceptance probability function input')


    def probability_acceptance(self, state_1: PolyhedralSpace, state_2: PolyhedralSpace) -> float:
        """
        This function compares two design states and calculates the product probability of acceptance to archive a
        design stage per the MOSA implementation

        :param state_1: First state (ie before grammar applied to design) being considered
        :param state_2: Second state (ie after grammar successfully applied to design) being considered
        :return: Probability acceptance to use in Metropolis-Hastings
        """
        # Find the objective function values of both design stages:
        old_objectives = self.calculate_objective_function_values(design_space=state_1)
        new_objectives = self.calculate_objective_function_values(design_space=state_2)

        # Calculate probabilities based on selected schedule:
        tot_prob = 0
        for old, new, temperature in zip(old_objectives, new_objectives, list(self.obj_temp_dict.keys())):
            # Calculate each probability for all objs.
            prob = self.p_accept(new=new, old=old, temperature=self.obj_temp_dict[temperature])
            # probs.append(prob)
            tot_prob += prob
        # product = reduce(mul, probs)
        final_probability = round(np.exp(tot_prob), self.numDecimals)
        return final_probability


    def metropolis(self, p_accept: float) -> bool:
        """
        Metropolis-Hastings algorithm implementation

        :param p_accept: Acceptance probability calculated
        :return: True: replace the design, False: Use probability test
        """
        ### PER THE MOSA ALGORITHM: Sometimes, p_accept will be greater than unity (1) for a given function. In this
        #                           case, we will always accept it. This occurs when the difference in the new and
        #                           old objective function values is very small (and we therefore want to explore this
        #                           region of the pareto, so we always accept these)
        if p_accept >= 1:
            return True

        # Otherwise, in the case of MINIMIZATION (which all of my objective functions are):
        self.p_accept_list.append(p_accept)  # Append to list for plotting later to verify algorithm
        if random() < p_accept:
            # If the random generated point is less than our acceptance probability, we return True to accept the state
            # Over time, the p_accept should go lower and lower meaning we accept less and less worse states.
            return True
        else:
            # Otherwise we do not accept the state and return False
            return False


    def store_objectives_worse_move(self, design_space: PolyhedralSpace) -> None:
        """
        Stores the objective function values when a worse move was selected (and we track here)

        design_space: Design space that is currently being evaluated
        """
        objective_vals = self.calculate_objective_function_values(design_space=design_space)
        for i in range(len(objective_vals)):
            self.tracked_objectives_at_Ti[i].append(objective_vals[i])
            self.tracked_objectives_all_temp[i].append(objective_vals[i])


    def quench_temperatures(self) -> None:
        """
        This function is in charge of cooling the temperature in MOSA with a few different cooling schedules as a
        potential use
        """
        if self.cooling_schedule == 'geometric' or self.cooling_schedule == 'Geoemtric':
            for i in range(len(self.objective_functions)):
                curTemp = self.obj_temp_dict[i]
                curTemp *= self.cooling_rate_geometric
                self.obj_temp_dict[i] = curTemp

        elif self.cooling_schedule == 'triki' or self.cooling_schedule == 'Triki':
            # Using a worst-case-scenario where all measured objectives are the same value, use a constant alpha of .5
            for i in range(len(self.objective_functions)):
                cur_objectives_list = self.tracked_objectives_at_Ti[i]
                if all(x == cur_objectives_list[0] for x in cur_objectives_list):
                    alpha = 0.5
                else:
                    alpha = (1 - ((self.obj_temp_dict[i] * self.delta_T) / np.std(cur_objectives_list) ** 2))
                    if alpha < 0.5:
                        alpha = 0.5  # We set a floor to the lowest the cooling rate can be
                # After calculating the new Triki cooling rate for this objective, we update tepmerature:
                curTemp = self.obj_temp_dict[i]
                curTemp *= alpha
                self.obj_temp_dict[i] = curTemp
                # Also, we reset the tracked objective list:
                self.tracked_objectives_at_Ti[i] = []  # Should reset the value of self.fX_obj_at_Ti  to []

        elif self.cooling_schedule == 'Huang' or self.cooling_schedule == 'huang':
            for i in range(len(self.objective_functions)):
                cur_objectives_list = self.tracked_objectives_at_Ti[i]
                if all(x == cur_objectives_list[0] for x in cur_objectives_list):
                    alpha = 0.5
                else:
                    alpha = max((0.5, np.exp(-0.7 * self.obj_temp_dict[i] / np.std(cur_objectives_list))))
                # After calculating the new Triki cooling rate for this objective, we update tepmerature:
                curTemp = self.obj_temp_dict[i]
                curTemp *= alpha
                self.obj_temp_dict[i] = curTemp
                # Also, we reset the tracked objective list:
                self.tracked_objectives_at_Ti[i] = []
        else:
            raise Exception('Invalid cooling schedule selected, options are: geometric, triki, huang')

        # Update temperature lists for plotting:
        for obj, temp in self.obj_temp_dict.items():
            curList = self.temp_tracker[obj]
            curList.append(temp)


    def store_pareto_for_animation(self, epoch: int) -> None:
        """
        This function takes in the current epoch and stores the required data for plotting an animation or video
        to see the growth / exploration of the design space from a pareto POV

        :param epoch: Current optimization epoch number to store the pareto data to
        """
        curData = list(self.MOSA_archive.keys())
        self.pareto_animation[epoch] = curData


    def begin_MOSA(self) -> None:
        """
        This method starts the optimization process and takes no input / has no output and simply writes out a file
        """
        # Initialize dictionaries for tracking and generate the initial shape for the input design space
        self.initialize_time_tracker_dict()
        self.design_space.generate_start_shape()

        ### STEP 0: (added this later, oops) Validate that the generated start shape actually passed all design
        # constraints to ensure that the optimization process can begin:
        self.design_constraints.create_constraint_dictionary(design_space=self.design_space)
        passed_all_constraints, constraint_times = self.design_constraints.check_constraints(
            design_space=self.design_space, constraint_set=self.constraint_set)
        if not passed_all_constraints:
            raise Exception('Invalid input design: The triangulated initial shape did not pass all design constraints'
                            ' of the problem, please review the constraints and initial triangulation of the problem'
                            ' definition. Please reach out to dev team if this is unclear.')

        ### STEP 1:  Initialize temperatures with all available grammars to get a better understanding of the
        #            design space w/ our grammars:
        if self.print_progress:
            print('Initializing temperatures for start state')
        self.initialize_temperatures()

        ### STEP 2 of MOSA: Use intelligent Return to Base to reset the reset the current design space along the pareto
        if self.print_progress:
            print('Performing initial return to base')
        self.return_to_base()

        ### STEP 3 of MOSA: Begin loops of Simulated Annealing with temperature full considered above
        inner_count = 0  # Used to track the internal loop count number so we know when to quench
        accepted_count = 0  # Used to track the number of accepted moves for quenching
        return_to_base_count = 0  # Used to track number of loops so we know when to return to base
        curEpoch = 1

        ### BEGIN MOSA (Currently a simple exit condition of Sim. Anneal)
        start_time = time.time()  # We use this to ensure a run doesn't go too long.
        while self.keep_annealing(start_time=start_time, numEpochs=curEpoch):
            # First we update the ramp extension values.
            extension_map, rotation_map = self.get_extend_and_rotation_values(epoch=curEpoch)
            # Step 1: Perturb the current design state into a sample state:
            new_grammar = self.grammar_set.pick_random_grammar()
            if 'Rotation' in new_grammar or 'Rotate' in new_grammar or 'rotation' in new_grammar:
                ext = rotation_map[new_grammar]
            else:
                ext = extension_map[new_grammar]

            # Create two copies of the design_space for comparison:
            active_design_state, candidate_design_state = deepcopy(self.design_space), deepcopy(self.design_space)

            # Apply the grammar to the candidate space and check the design constraints:
            grammar_applied_successfully, grammar_time = self.grammar_set.call_grammar_function(
                grammar_selected=new_grammar,
                design_space=candidate_design_state,
                extension_value=ext)
            if grammar_applied_successfully:
                # Update the design constraint values and calculate if design is compliant:
                self.design_constraints.update_params(design_space=candidate_design_state)
                passed_all_constraints, constraint_times = self.design_constraints.check_constraints(
                    design_space=candidate_design_state, constraint_set=self.constraint_set,
                    compare_space=self.design_space)
            else:
                passed_all_constraints = False  # If grammar wasn't applied don't both calculating these constraints

            # For tracking we update the time_tracker for the given grammars:
            if grammar_time > self.time_tracker_dict[new_grammar]:
                self.time_tracker_dict[new_grammar] = grammar_time
            self.update_constraint_time_tracker(new_times=constraint_times)  # Update tracking of times...

            if passed_all_constraints and grammar_applied_successfully:
                # If we pass the constraints check AND the grammar applied successfully, we attempt archiving:
                archive_design = self.attempt_to_archive(candidate_design=candidate_design_state)

                if archive_design:
                    # If selected for archiving, the candidate_design_state is set to our new design state:
                    del active_design_state
                    self.update_design_space(new_design_space=candidate_design_state)
                else:
                    # Otherwise, we use a probability test to determine acceptance:
                    p_accept = self.probability_acceptance(state_1=active_design_state, state_2=candidate_design_state)
                    accept_new_state = self.metropolis(p_accept=p_accept)
                    """
                    Note: Overall, this selection criterion used is a weak point of this MOSA and it is only 
                          recommended for ~two/three objectives. Further research should be conducted here to test the 
                          generalizable to higher number of objective functions. I won't be going higher than 3 tho as 
                          the actual Pareto is obfuscated due to visualization in 4D or higher being tricky to understand
                    """
                    if accept_new_state:
                        # If we accept via Metropolis then we set the new state to the active state and store objectives
                        del active_design_state
                        self.update_design_space(new_design_space=candidate_design_state)
                        self.store_objectives_worse_move(design_space=candidate_design_state)
                        self.accepted_via_probability += 1
                        self.total_checked += 1
                        accepted_count += 1  # If we probabilistically accept, we iterate accepted_count
                    else:
                        # If we do not accept via Metropolis, we just continue along with X_N1 as the active state
                        self.total_checked += 1

                    if self.total_checked > 1:  # We only start counting once we have passed a constraint check
                        self.acceptance_tracker.append(self.accepted_via_probability / self.total_checked)
            inner_count += 1  # Always increment inner_count even if constraints were not passed

            ### PERIODICALLY: Return to base strategy to explore Pareto
            return_to_base_count += 1  # Increment the return to base count to determine when we should next return:
            if return_to_base_count == self.N_Bi:
                # If we are returning to base, then we call the new state and set it to XN1:
                self.return_to_base()
                return_to_base_count = 0  # Reset this count here

            ### PERIODICALLY: Annealing of Temperature per van Laarhoven and Aarts
            if inner_count >= self.NT2 or accepted_count >= self.Na:
                if inner_count >= self.NT2 and self.print_progress:
                    print(f'inner count of {inner_count} exceeded NT2 of {self.NT2} at epoch: {curEpoch}')
                    print(f'The value of accepted count was {accepted_count} at epoch: {curEpoch}')
                elif accepted_count >= self.Na and self.print_progress:
                    print(f'accepted count of {accepted_count} exceeded Na of {self.Na} at epoch: {curEpoch}')
                    print(f'The value of inner_count was {inner_count} at epoch: {curEpoch}')
                self.quench_temperatures()
                # Now we update the acceptance_per_epoch tracker:
                self.acceptance_per_epoch.append(np.mean(self.acceptance_tracker))
                self.acceptance_tracker = []  # Reset to an empty list

                # After quenching we set the counters for BOTH inner_count and accepted_count to 0 for the next temperature
                inner_count, accepted_count = 0, 0
                self.store_pareto_for_animation(epoch=curEpoch)
                curEpoch += 1  # When quenching, we add 1 to epoch.


        # After the while loop concludes / exits, we return the ARCHIVED dataset so we can access any of these values
        end_time = time.time()
        self.sim_time_seconds = end_time - start_time

        # Finally, we save the output by calling the function:
        self.save_output()


    def save_output(self) -> None:
        """
        This function is called at the end of a MOSA operation to save the output file

        :return: Nothing, will simply write a dill binary file containing class information at final iteration
        """
        # I am using the 'aj1' extension because I want to track which versions of my software can generate which
        # types of solutions and plots. For example, if I add in features later on I want to differentiate those types
        # of simulation output.
        def is_valid_filename(filename):
            import re
            # Define a regular expression for a valid filename (letters, numbers, and underscores only)
            pattern = re.compile(r'^[a-zA-Z0-9_]+$')

            # Check if the filename matches the pattern
            if pattern.match(filename):
                return True
            else:
                return False

        # Check if the save path exists, if not make it:
        if not os.path.exists(self.SAVE_PATH):
            os.mkdir(self.SAVE_PATH)

        # Create file:
        if is_valid_filename(filename=self.SAVE_NAME_NO_EXTENSION):
            save_name = self.SAVE_NAME_NO_EXTENSION + '.aj1'
        else:
            warnings.warn("Detected an invalid filename, therefore resetting to randomized filename.")
            save_name = 'MOSA_Export_Data.aj1'

        output_file = os.path.join(self.SAVE_PATH, save_name)
        copy_of_data = deepcopy(self)
        with open(output_file, 'wb') as f:
            dill.dump(copy_of_data, f)

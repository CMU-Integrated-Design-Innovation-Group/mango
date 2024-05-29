"""
A.J. Vetturini
IDIG and MMBL
Carnegie Mellon University

Currently this optimizer is only available for the PolyhedralSpace; future work will potentially need to incorporate
different spaces and update this class to allow different spaces.

This script incorporates the shape annealing algorithm as developed in Cagan, J., and W. J. Mitchell.
“Optimally Directed Shape Generation by Shape Annealing.” Environment and Planning B: Planning and Design 20,
no. 1 (1993): 5–12.
"""

# Import Modules
from dataclasses import dataclass, field
from copy import deepcopy
from typing import List, Union
import warnings
from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
from mango.grammars.origami_grammars import GrammarSet, CustomGrammarSet
from mango.optimization_features import design_constraints, objective_function
from mango.utils.mango_math import *
import time
from random import random, seed, randint
import os
import dill

@dataclass
class ShapeAnneal(object):
    """
    The optimizer is the "heart and soul" of the framework. It contains the logic for applying grammars and controlling
    the optimization-driven generative process. There are various hyperparameters that can be tuned which establish
    the depth of search and time spent searching.

    This data class contains the Shape Annealing algorithm implementation and is skewed towards a single objective
    optimization task. It can be used with a weighted sums / epsilon constraint approach to multiobjective, but it
    is not necessarily an efficient way as compared to the multiobjective shape annealing implementation.


    Parameters
    ------------
    design_space : Initial PolyhedralSpace representing the design which is to be optimized
    grammars : A singular GrammarSet or a list of GrammarSets to use in the generative process
    design_constraints : The PolyhedralDesignConstraint object which constraint the optimizer
    objective_function : Objective function to be minimized
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

    Simulated Annealing Hyperparameters that can be controlled
    Note that the default parameters here won't necessarily give a "good" result, they must be tuned!
    ------------
    n : Number of inner loop mutations (default: 50)
    limit : Set limit / # of accepted shapes within inner most loop of SA (default: 15)
    n_stop : Max number of iterations to prevent an infinite loop (default: 100000)
    T_min : Minimal temperature to end the annealing process (default: 1e-8)
    max_time_of_optimization_minutes : Max time spent generating a design (default: 60 minutes)
    max_number_of_epochs : Max number of epochs to run the generative study for (default: 50)
    random_seed : Specified random seed (default: Random value)
    random_walk_steps : Number of random walk steps to initialize the temperature of the space (default: 2000)
    cooling_schedule : Cooling schedule used during simulated annealing (default: huang / HRSV)
    cooling_rate_geometric : Cooling rate to use if using geometric cooling_schedule (default: 0.9)
    delta_T : Triki annealing schedule adaptation rate, only used if using triki cooling_schedule (default: 0.8)

    Optimization Convergence Parameters:
    ------------
    FYI:  Because this work was published in a bio-focused field, we do NOT use a convergence set by default, but rather
          let the optimization process carry on for a set time. This is because this framework is really intended for
          design conceptualization, and fine tuning a design is not going to be accomplished with this framework.
    use_convergence : Boolean to use convergence or whether to just run for a set amount of time (default : FALSE)
    minimal_convergence_size : Integer number of # of designs to consider when calculating convergence (default : 25)
    convergence_threshold : Floating point threshold for determining design convergence (default = 0.01)
    """
    # Input parameters that must be passed in:
    design_space: PolyhedralSpace
    grammars: Union[GrammarSet, list]  # Pass a list of any grammar sets used.
    design_constraints: design_constraints
    objective_function: objective_function.ObjectiveFunction
    SAVE_PATH: str  # This is where the output dill file is saved to
    SAVE_NAME_NO_EXTENSION: str
    print_progress: bool = True  # Print out messages to monitor the optimization process
    # Can be changed to use different constraint sets within a design_constraints obj, just input names of constraint
    # functions found within design_constraints to be called
    constraint_set: List[str] = field(default_factory=lambda: ['Outside Design Space', 'Vertex in Excluded',
                                                               'Edge in Excluded', 'Invalid Edge Length',
                                                               'Invalid Scaffold Length', 'Invalid Face Angle',
                                                               'Broken preserved edge', 'Intersecting edges',
                                                               'Intersecting faces'])

    # Simulated Annealing hyperparameters user can control:
    n: int = 50  # Number of inner loop mutations
    limit: int = 15  # Set limit / # of accepted shapes within inner most loop of SA
    n_stop: int = 100000  # Max number of iterations to prevent an infinite loop
    T_min: float = 1e-8  # Minimal temperature to end the annealing process
    max_time_of_optimization_minutes: int = 60
    max_number_of_epochs: int = 50
    random_seed: int = None  # Random seed a user can specify

    # Various input parameters that can be changed:
    extension_value_default: float = 0.34  # Default value if ramp is not used
    extension_ramp: dict = field(default_factory=dict)
    rotation_value_degrees_default: float = 1
    rotation_ramp: dict = field(default_factory=dict)
    numDecimals: int = 5

    # Parameters use can control to determine convergence
    use_convergence: bool = False
    minimal_convergence_size: int = 25
    convergence_threshold: float = 0.01

    ## Hyperparameters that can be tuned by user:
    random_walk_steps: int = 2000  # Number of random walk steps to initialize the temperature of the space
    cooling_schedule: str = 'Huang'  # Cooling schedule used defined by user

    # Cooling Schedule Parameters depending on schedule used (Huang is recommended)
    cooling_rate_geometric: float = 0.8  # The cooling schedule for geometric (only used if Geoemtric schedule is used)
    delta_T: float = 0.8  # Used with the Triki annealing schedule and controls adaptation. May need to modify this.


    def __post_init__(self):
        self.sim_time_seconds = 0.0
        self.final_time = None
        self.time_tracker_dict = {}
        self.design_evolution = {}
        self.random_walk_objective_values = []
        self.temperature_values = []
        self.objectives_during_annealing = []
        self.optimization_objective_values = []
        self.tracked_objectives_at_Ti = []
        self.T = 0  # initial temperature
        self.scaffold_needed = []

        if self.random_seed is None:
            self.random_seed = randint(1, 10000000)
            print(f'No random seed specified, using the following seed number: {self.random_seed}')
        seed(self.random_seed)  # Set the random seed
        np.random.seed(self.random_seed)

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

        self.verify_input_values()
        # Validate ramp conditions:
        if self.extension_ramp != {}:
            # If the user specifies an extension_ramp, then we must validate the values to ensure it is a proper
            # data structure:
            for grammar, ramp in self.extension_ramp.items():
                if grammar not in self.grammar_set.grammar_names:
                    raise Exception(f'Invalid ramp specified for {grammar} as this grammar is not in the list of grammars '
                                    f'specified to the problem.')
            # If we validate all ramps w/o error, then we just set the flag to true:
            self.use_extension_ramp = True
        else:
            self.use_extension_ramp = False

        # Repeat above for the rotation_ramp.
        if self.rotation_ramp != {}:
            for grammar, ramp in self.rotation_ramp.items():
                if grammar not in self.grammar_set.grammar_names:
                    raise Exception(f'Invalid ramp specified for {grammar} as this grammar is not in the list of grammars '
                                    f'specified to the problem.')
            self.use_rotation_ramp = True
        else:
            self.use_rotation_ramp = False


    def verify_input_values(self):
        """ Function that simply validates the input parameters before starting an optimization process. """
        if self.cooling_schedule not in ['geometric', 'Geometric', 'triki', 'Triki', 'Huang', 'huang', 'HRSV']:
            raise Exception('Invalid cooling schedule input, only values are: geometric, triki, huang')


    def initialize_time_tracker_dict(self):
        """ Initialize a dictionary to track time for internal development purposes """
        # first create keys for the grammars and constraints. We use "0" since a time-to-execute will never be less than
        # 0 and these are just tracking the execution-times of these functions.
        for grammar in self.grammar_set.grammar_names:
            self.time_tracker_dict[grammar] = 0

        for constraint in self.design_constraints.names:
            self.time_tracker_dict[constraint] = 0

        # Add in some defaults for the objective function evaluation timer
        nameMax = self.objective_function.name + "_max_time"
        nameMin = self.objective_function.name + "_min_time"
        # Initialize the max and min values to values that will get over-written during the tracking process
        self.time_tracker_dict[nameMax] = -1
        self.time_tracker_dict[nameMin] = 1e6


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


    def update_design_space(self, new_design_space: PolyhedralSpace) -> None:
        """ This method updates the active design_space object whenever it needs to be completely updated """
        self.design_space = new_design_space


    def calculate_objective_function_value(self, design_space: PolyhedralSpace) -> float:
        """
        This function is responsible for evaluating the design_space using the input objective function. It also records
        the execution times for internal development purposes.

        :param design_space: Design space that is currently being evaluated
        :returns: Floating point value of objective function value
        """
        all_edge_lengths = calculate_design_edge_lengths(graph=design_space.design_graph)
        input_params = design_space.calculate_input_parameters(edge_lengths=all_edge_lengths,
                                                               routing_algorithm=self.design_constraints.scaffold_routing_algorithm)
        nameMax = self.objective_function.name + "_max_time"
        nameMin = self.objective_function.name + "_min_time"
        st_t = time.time()
        obj_val = self.objective_function.evaluate_function(input_params)
        end_t = time.time()
        elapsed = end_t - st_t
        # RECORD TIMES
        if elapsed > self.time_tracker_dict[nameMax]:
            self.time_tracker_dict[nameMax] = elapsed
        if elapsed < self.time_tracker_dict[nameMin]:
            self.time_tracker_dict[nameMin] = elapsed
        # RECORD OBJECTIVE FUNCTION VALUE
        return obj_val


    def calculate_and_store_objective_function_values(self, design_space: PolyhedralSpace, random_walk: bool) -> None:
        """
        This function is called so that the objective function is stored to analyze the optimization performance.

        :param design_space: Design space that is currently being evaluated
        :param random_walk: Boolean determining if this is a random_walk step or not (False)
        """
        obj_value = self.calculate_objective_function_value(design_space=design_space)
        if random_walk:
            self.random_walk_objective_values.append(obj_value)
        else:
            self.optimization_objective_values.append(obj_value)
            self.scaffold_needed.append(design_space.input_values['estimated_scaffold_length'])


    def store_design_for_animation(self, epoch: int) -> None:
        """
        Method to simply store a design at every epoch in the simulated annealing process to visualize how a design
        evolved during the generative process.

        :param epoch: Current epoch of the optimizer to assign the stored design to
        """
        # First evaluate the objective:
        obj_value = self.calculate_objective_function_value(design_space=self.design_space)
        # Store in the dictionary in the format EPOCH : (objective_valuation, design_space object)
        self.design_evolution[epoch] = (obj_value, self.design_space)


    def initialize_temperature(self) -> None:
        """
        This function performs a random walk through the objective space to automatically assign the temperature used
        in the acceptance criterion for the first epoch.
        """
        # In a random walk we will ALWAYS accept a move, but we need a counter and while loop since we can still
        # violate constraints
        accepted_moves = 0
        while_loop_counter = 0  # Prevent infinite loop

        # Get the extension value from the ramp. This will return the defaults if a ramp is not used as a note!
        extension_map, rotation_map = self.get_extend_and_rotation_values(epoch=0)

        initial_design_to_return = deepcopy(self.design_space)
        # Random Walk:
        while accepted_moves < self.random_walk_steps:
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
                    design_space=candidate_design_state, constraint_set=self.constraint_set, compare_space=self.design_space)
            else:
                # If the grammar didn't apply just set this to False. Doesn't need to be this way but ;-)
                passed_all_constraints = False

            if grammar_applied_successfully and passed_all_constraints:
                # If the grammar is applied successfully and the constraints are all passed then this new state replaces
                # the current design_space
                self.calculate_and_store_objective_function_values(design_space=candidate_design_state,
                                                                   random_walk=True)
                self.update_design_space(new_design_space=candidate_design_state)
                accepted_moves += 1

            #  No matter what, we delete this candidate from memory so it re-instantiates in while loop for peace of
            #  mind regarding objects. If the design_space was accepted in if above, then this all works out!
            del candidate_design_state

            # For the sake of no infinite loops, I use a while loop counter here to ensure it doesn't get stuck:
            while_loop_counter += 1
            if while_loop_counter > (10 * self.random_walk_steps):  # Prevent infinite loop
                warnings.warn(f'Unable to apply {self.random_walk_steps} without violating design constraints which may '
                              f'lead to an inferior solution set, be aware!')
                accepted_moves = self.random_walk_steps + 1  # Just break the loop condition...

        # After the while loop has ended, we set the value of temperature to the std. dev of the objectives tracked:
        newTemp = round(float(np.std(self.random_walk_objective_values)), self.numDecimals)
        if np.isclose(newTemp, 0):
            raise Exception(f'Could not find any variation in the objective function valuation for the'
                            f'specified function {self.objective_function.name}. Please increase random_walk_steps or '
                            f'verify the objective function is evaluating properly!')
        self.T = newTemp  # Assign the std. dev to the initial temperature value
        self.temperature_values.append(self.T)  # Store this initial value in the list
        # Return to the initial state to start the optimization process from
        self.update_design_space(new_design_space=initial_design_to_return)
        # Store the initial design now that temperature is initialized:
        self.store_design_for_animation(epoch=0)


    def keep_annealing(self, start_time: float, cur_epoch: int) -> bool:
        """
        This function determines when the SA algorithm should end based on a length of time or if the temperature
        hyperparameter gets too low

        :param cur_epoch: Current epoch of SA
        :param start_time: Start time of how long the inner while loop has been running
        :return: True = Keep running, False = END
        """
        elapsed_time = time.time() - start_time
        hours, remainder = divmod(elapsed_time, 3600)
        minutes, seconds = divmod(remainder, 60)
        time_formatted = "{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds)

        if self.use_convergence:
            # If the user signals to use a convergence criterion we use this to determine if we should end
            ## First grab the last N (N=minimal convergence size) datapoints:
            if len(self.objectives_during_annealing) < self.minimal_convergence_size:
                # If we do not have enough datapoints, we can NOT converge
                return True
            obj_vals = self.objectives_during_annealing[-self.minimal_convergence_size:]
            changes = [abs(obj_vals[i] - obj_vals[i - 1]) for i in range(1, len(obj_vals))]
            if all(change < self.convergence_threshold for change in changes):
                # If the convergence is beneath the threshold for all changes, we return False signalling "end process"
                return False

            # However we still check max number of epochs to prevent an infinite loop:
            if cur_epoch == self.max_number_of_epochs:
                # If the current epoch is the value of the max_number, then we end
                if self.print_progress:
                    print(f'Shape Annealing ended due to the max number of epochs being met: '
                          f'{self.max_number_of_epochs}')
                self.final_time = time_formatted
                return False

            # Otherwise return True to signal "keep optimizing"
            return True


        else:
            ### First check if temperatures are low enough to end the run:
            if self.T < self.T_min:
                if self.print_progress:
                    print(f'Shape Annealing ended due to temperatures dropping below threshold of {self.T_min}')
                self.final_time = time_formatted
                return False

            ### Next check if max time has elapsed
            if (time.time() - start_time) >= (self.max_time_of_optimization_minutes * 60):
                # If the current elapsed time is greater than the threshold in seconds (why i multiply by 60 above), we end
                if self.print_progress:
                    print(f'Shape Annealing ended due to max time input of {self.max_time_of_optimization_minutes} minutes')
                self.final_time = time_formatted
                return False

            if cur_epoch == self.max_number_of_epochs:
                # If the current epoch is the value of the max_number, then we end
                if self.print_progress:
                    print(f'Shape Annealing ended due to the max number of epochs being met: '
                          f'{self.max_number_of_epochs}')
                self.final_time = time_formatted
                return False

            else:
                # Otherwise we return True meaning we keep looping
                return True


    def metropolis(self, active_design: PolyhedralSpace, candidate_design: PolyhedralSpace) -> bool:
        """
        This function determines if a candidate_design should be accepted and replace the active design. The
        metropolis algorithm is used to determine this acceptance:

        :param active_design: S1 state in Sim. Annealing
        :param candidate_design: S2 state in Sim. Annealing
        :return: True: replace the design, False: Use probability test
        """
        s1 = self.calculate_objective_function_value(design_space=active_design)
        s2 = self.calculate_objective_function_value(design_space=candidate_design)
        # In Metropolis-Hastings if the new is less that the previous we accept unconditionally for minimization
        if s2 < s1:
            self.tracked_objectives_at_Ti.append(s2)  # Store the objective value for annealing schedule if accepted
            self.objectives_during_annealing.append(s2)  # Store objective also for all objective func. values tracking
            return True
        else:
            # Otherwise, we use a probability acceptance and a random value to determine if we should accept:
            test_value = np.exp(-abs(s2 - s1) / self.T)
            if random() < test_value:
                # If the random value is less than this test value, we accept the worse move
                self.tracked_objectives_at_Ti.append(s2)  # Store the objective value for annealing schedule if accepted
                self.objectives_during_annealing.append(s2)
                return True
            else:
                # Otherwise, we do not accept the candidate design
                return False


    def quench_temperatures(self, cur_epoch: int) -> None:
        """
        This function is in charge of cooling the temperature in SA with a few different cooling schedules as a
        potential use.

        WIP: Should add some simple regex, this is a dumb way of checking strings used.
        """
        # First we store the design for the animation at the current epoch:
        self.store_design_for_animation(epoch=len(self.objectives_during_annealing))  # This aligns chart better?

        if self.cooling_schedule == 'geometric' or self.cooling_schedule == 'Geoemtric':
            self.T *= self.cooling_rate_geometric

        elif self.cooling_schedule == 'triki' or self.cooling_schedule == 'Triki':
            # Using a worst-case-scenario where all measured objectives are the same value, use a constant alpha of .5
            if all(x == self.tracked_objectives_at_Ti[0] for x in self.tracked_objectives_at_Ti):
                alpha = 0.5
            else:
                alpha = (1 - ((self.T * self.delta_T) / np.std(self.tracked_objectives_at_Ti) ** 2))
                if alpha < 0.5:
                    alpha = 0.5  # We set a floor to the lowest the cooling rate can be
            # After calculating the new Triki cooling rate for this objective, we update temperature:
            self.T *= alpha
            # Also, we reset the tracked objective list:
            self.tracked_objectives_at_Ti = []  # Should reset the value of self.fX_obj_at_Ti  to []

        elif self.cooling_schedule == 'Huang' or self.cooling_schedule == 'huang':
            if all(x == self.tracked_objectives_at_Ti[0] for x in self.tracked_objectives_at_Ti):
                alpha = 0.5
            else:
                alpha = max((0.5, np.exp(-0.7 * self.T / np.std(self.tracked_objectives_at_Ti))))
            # After calculating the new Triki cooling rate for this objective, we update temperature:
            self.T *= alpha
        else:
            raise Exception('Invalid cooling schedule selected, options are: geometric, triki, huang')

        self.temperature_values.append(self.T)
        self.tracked_objectives_at_Ti = []


    def begin_annealing(self) -> None:
        """
        This is the top-level method that a user calls to begin the generative process. Nothing is passed in or returned
        """
        # Generate start shape and initialize time tracking dictionaries for logging optimizer performance
        self.initialize_time_tracker_dict()
        self.design_space.generate_start_shape()

        ### STEP 0: (added this later, oops) Validate that the generated start shape actually passed all design
        # constraints to ensure that the optimization process can begin. Here, we also intialize the design constraint
        # tracking for computational efficiency
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
        self.initialize_temperature()

        start_time = time.time()
        curEpoch = 1
        while self.keep_annealing(start_time=start_time, cur_epoch=curEpoch):
            success = 0  # Count of successive state applications

            # Inner loop of n steps to modify design:
            for _ in range(self.n):
                active_design_state, candidate_design_state = deepcopy(self.design_space), deepcopy(self.design_space)
                extension_map, rotation_map = self.get_extend_and_rotation_values(epoch=curEpoch)
                # Randomly pick and apply a grammar:
                new_grammar = self.grammar_set.pick_random_grammar()
                if 'Rotation' in new_grammar or 'Rotate' in new_grammar or 'rotation' in new_grammar:
                    ext = rotation_map[new_grammar]
                else:
                    ext = extension_map[new_grammar]

                grammar_applied_successfully, grammar_time = self.grammar_set.call_grammar_function(
                    grammar_selected=new_grammar,
                    design_space=candidate_design_state,
                    extension_value=ext)  # Only pass in the specified value
                if grammar_applied_successfully:
                    # Update the design constraint values and calculate if design is compliant:
                    self.design_constraints.update_params(design_space=candidate_design_state)
                    passed_all_constraints, constraint_times = self.design_constraints.check_constraints(
                        design_space=candidate_design_state, constraint_set=self.constraint_set, compare_space=self.design_space)
                else:
                    passed_all_constraints = False  # If grammar wasn't applied don't both calculating these constraints

                # For tracking we update the time_tracker for the given grammars:
                if grammar_time > self.time_tracker_dict[new_grammar]:
                    self.time_tracker_dict[new_grammar] = grammar_time

                # Now check if we passed all constraints and applied the grammar successfully to check to see if we
                # should accept this design state:
                if passed_all_constraints and grammar_applied_successfully:
                    # If this condition is met we evaluate both the active and candidate designs:
                    accept_design = self.metropolis(active_design=active_design_state,
                                                    candidate_design=candidate_design_state)
                    if accept_design:
                        # If the design is accepted then we update the design state and calculate and store the
                        # objective function value
                        del active_design_state
                        self.update_design_space(new_design_space=candidate_design_state)
                        self.calculate_and_store_objective_function_values(design_space=candidate_design_state,
                                                                           random_walk=False)
                        success += 1

                # Now we check to see this inner loop should break early:
                if success > self.limit:
                    if self.print_progress:
                        print(f'Quenching temperature due to successes {success} exceeding limit {self.limit}')
                    break

            # After the for loop we check if successes is 0. If so then the equilibrium value is met and we end!
            # Note: This is a bad convergence criterion with the use of a ramp element / in general imo
            #       I'd rather use the L1 or just for a # of epochs, this isn't inherently always true.
            """if success == 0:
                break"""

            # Now just quench temperatures and increment the epoch counter
            self.quench_temperatures(cur_epoch=curEpoch)
            curEpoch += 1

        # After the while loop ends we export the design:
        end_time = time.time()
        self.sim_time_seconds = end_time - start_time

        # Finally, we save the output by calling the function:
        self.save_output()


    def save_output(self) -> None:
        """
        This function is called at the end of a SA algorithm to save the output file

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
            os.makedirs(self.SAVE_PATH)

        # Create file:
        if is_valid_filename(filename=self.SAVE_NAME_NO_EXTENSION):
            save_name = self.SAVE_NAME_NO_EXTENSION + '.aj1'
        else:
            warnings.warn("Detected an invalid filename, therefore resetting to randomized filename.")
            save_name = 'Single_Objective_Export_Data.aj1'


        output_file = os.path.join(self.SAVE_PATH, save_name)
        copy_of_data = deepcopy(self)
        del copy_of_data.design_constraints.names_to_functions  # This messes up this dill object
        with open(output_file, 'wb') as f:
            try:
                dill.dump(copy_of_data, f)
            except TypeError as e:
                # Check if the error is due to a specific type and skip over it
                print(f'Error exporting file: {e}')
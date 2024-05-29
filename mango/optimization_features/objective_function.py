"""
A.J. Vetturini
IDIG and MMBL
Carnegie Mellon University

This script simply contains the dataclass for creating an ObjectiveFunction to use to guide the search process for
a generative study.
"""
from dataclasses import dataclass
from typing import Callable

@dataclass
class ObjectiveFunction(object):
    """
    This class lets the user define a custom function which can be assigned as the objective of the generative process.
    This is a slight WIP as I may want to convert this to a Python wrapper, but I am not sure about the logistics (or if
    this is even a good idea). This class is very similar to the CustomDesignConstraint class.

    The function can really be defined as anything, however it should be noted that there are no checks here to see how
    "good" a constraint is or is not, and that is at the discretion of the user.

    Parameters
    ------------
    name : String unique identifies of the constraint name
    objective_equation : A Python function defining the objective. For example, if you have a constraint as:
                        def bar(input_vars, extra_params):
                            objective_function_code

                        Then you would pass in ObjectiveFunction(name='bar', design_constraint=bar,
                                                                 extra_params=extra_params)
    extra_params : A dictionary of constants / values that a user may want to use in the custom constraint
    numDecimals : Floating point value to control precision. Allowing too high of precision, depending on the sampling
                  sampling, can lead to the optimizers getting stuck in local minima around a certain point.
    """
    name: str  # Name of the objective function to show on plots
    objective_equation: Callable = None  # This is the actual equation (which can ONLY contain ONE term) in LaTeX formatting for sympy
    extra_params: dict = None  # Extra values that the user may pass in
    # Constant value that a user might want to use to properly "scale" the objective function. Note that each Objective
    # Function class will only have "one" component so there is only 1 constant value.
    numDecimals: int = 3  # The higher this value the more MOSA can get trapped in local minima

    def evaluate_function(self, input_parameters):
        # If we are using a "simple" (or pre-programmed) function, we have to set it as a callable function:
        try:
            return round(self.objective_equation(input_parameters, self.extra_params), self.numDecimals)
        except TypeError:
            return round(self.objective_equation(input_parameters), self.numDecimals)
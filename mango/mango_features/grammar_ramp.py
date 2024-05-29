"""
A.J. Vetturini
IDIG and MMBL
Carnegie Mellon University

This script contains the ability to create ramp elements for both the edge rotation and the length extension rules.
These ramps are specifically used for the length extension or edge rotation rules, future research may seek to optimize
these ramps automatically based on the scale of the input.
"""
from dataclasses import dataclass
import numpy as np
import warnings

@dataclass
class Ramp(object):
    """
    This Ramp class is capable of manipulating the extension / rotation value of the grammars during the generative
    process. Generally, earlier in the search we make larger design modifications and over time we switch over to a
    more refined search by making smaller design modifications down to a minmal value.

    Parameters
    ------------
    unique_id : Center of sphere defined by a numpy array in [X Y Z] (nm).
    max_value : floating point of the sphere (nm)
    min_value : Number of decimals to round positional values to
    max_number_epochs : Integer number defining how many epochs a study might run for to allocate step sizes
    min_steps_at_min : Integer number of a minimal # of steps to use the min_value for the ramp element
    precision : integer that is automatically specified to 2. You can lower the precision if desired but consider DNA!
    """
    unique_id: str  # A simple identifier / name for the ramp object
    max_value: float  # User should be careful with too large, we currently employ an error to prevent this behavior
    min_value: float  # E.g. 0.34nm since we only need precision down to the level
    max_number_epochs: int  # So we know how to ramp based on Epoch
    min_steps_at_min: int = 10  # Number of steps where the extend value is the lowest possible (0.34 nm)
    precision: int = 2  # Level of precision for the ramp, 2 is good for DNA since BDNA bprise is 0.34nm and
                        # we don't need more precision than that.

    def __post_init__(self):
        if self.max_value > (20 * self.min_value): # If desired, just raise the value of "20" but please note
                                                   # this in any publication supplmentary-type materials.
            raise Exception('Too large of a discrepancy between the max and min value for a ramp.')

        self.ramp = self.create_ramp()  # Make this function and things later.

    @staticmethod
    def distribute(ct: int, values: np.ndarray) -> list:
        """ This static method simply distributes the ramp elements given the parameters defined at class level """
        num_bins = len(values)
        ideal_value_per_bin = ct // num_bins
        remainder = ct % num_bins

        # Sort the values in ascending order
        sorted_values = sorted(values)

        distributed_values = []
        for _ in sorted_values:
            bin_value = ideal_value_per_bin
            if remainder > 0:
                bin_value += 1
                remainder -= 1
            distributed_values.append(bin_value)

        return distributed_values


    def create_ramp(self) -> np.ndarray:
        """
        Note: This is an area of interest for the future, here I use a simple ramp that worked via design iteration.
              Future studies can observe the effects of different ramps (and maybe influence repeatability and depth
              of search given the stochastic nature of the optimizers).

        This class creates the actual ramp array which can be referenced by the optimizer to determine the extension
        or rotation value during that design iteration.
        """
        max_mult = round(self.min_value * np.round(self.max_value / self.min_value), self.precision)
        if max_mult != self.max_value:
            warnings.warn(f'The maximal value for the {self.unique_id} ramp of {self.max_value} was not an even step '
                          f'size away from the min value of {self.min_value}, therefore using a value of {max_mult} '
                          f'instead.')
        extension_values = np.linspace(self.min_value, max_mult, num=int(np.round(self.max_value/self.min_value)))
        ramp = np.zeros(self.max_number_epochs)
        idx = len(extension_values)
        ct = self.max_number_epochs - self.min_steps_at_min  # Number of "steps" we have to assign to non-min vals
        if ct > 0:
            # If ct is larger than the number of value in ramp, then we try and assign as many "even" time at each of
            # values. However, we prefer to spend more time at the smaller values depending on the length discrepancy
            # between ct and extension values
            # This distribute logic bakes in if ct < len(extension_values) and we will only use the smallest values
            distributions = self.distribute(ct=ct, values=extension_values)
            flipped_distributions = distributions[::-1]

            ct = 0
            for i in flipped_distributions:
                for _ in range(i):
                    ramp[ct] = extension_values[idx-1]
                    ct += 1
                idx -= 1  # We lower idx to lessen the input value into flipped distributions

        # If ct is not bigger than zero we do the logic below. This logic still "works" to finish off the ramp as
        # edited above as in this case, all we need to do is replace all remaining "zeros" in ramp with the value
        # of min_value:
        mask = ramp == 0.0
        ramp[mask] = self.min_value

        return ramp


    def current_ramp_value(self, epoch: int) -> float:
        """
        This function returns the ramp value to use based on the current epoch #

        :param epoch: Whatever epoch # we are at
        :return: The ramp value to use in the optimizer
        """
        try:
            return self.ramp[epoch]
        except IndexError:
            # In case an issue, just use the last value (ie smallest value)
            return self.ramp[-1]

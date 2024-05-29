"""
A.J. Vetturini
IDIG and MMBL
Carnegie Mellon University

This script contains the DNA constant values (e.g. bp rise)
"""
from dataclasses import dataclass

@dataclass(frozen=True)
class BDNA(object):
    diameter: float = 2.25  # Presumed diameter of BDNA used in DNA origami
    pitch_per_rise: float = 0.34  # 0.34nm / basepair is the presumed rise for DNA origami
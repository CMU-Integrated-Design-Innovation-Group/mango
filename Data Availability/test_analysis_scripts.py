"""
This script valildates the optimizer_analysis functions
"""
from InternalFeatures import edge_cutoff_constraint, porosity_objective, cylinder_volume

from mango.visualizations.optimizer_analysis import ShapeAnneal_Analysis
single_objective_file = "./data_output/InternalFeatures2/10000_ExcludedFalse_Cutoff50_internal_study_1.aj1"

sa = ShapeAnneal_Analysis(results_path=single_objective_file)
sa.plot_objective_trace()
sa.show_fig()
sa.plot_temperature_trace()
sa.show_fig()

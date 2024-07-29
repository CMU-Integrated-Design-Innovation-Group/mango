"""
This script valildates the optimizer_analysis functions
"""
from mango.visualizations.optimizer_analysis import MOSA_Analysis, ShapeAnneal_Analysis
single_objective_file = "/Users/kodak/Desktop/mango/mango/tests/test_outputs/shape_annealing_test.aj1"
multi_objective_file = "/Users/kodak/Desktop/mango/mango/tests/test_outputs/MOSA_test.aj1"

sa = ShapeAnneal_Analysis(results_path=single_objective_file)
sa.plot_objective_trace()
sa.show_fig()
sa.plot_temperature_trace()
sa.show_fig()

mo = MOSA_Analysis(results_path=multi_objective_file)
mo.plot_all_objective_func_values()
mo.show_fig()
mo.plot_acceptance_probabilities()
mo.show_fig()
mo.plot_temperatures_over_time()
mo.show_fig()
mo.plot_final_2D_paretos()
mo.show_fig()
mo.create_pareto_animation_figure()
mo.show_fig()
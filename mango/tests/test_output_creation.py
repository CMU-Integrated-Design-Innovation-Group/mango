from mango.visualizations.mango_output_visualization import VisualizationObject
from mango.utils.design_io import mass_design_export

# Single objective output:
'''results_filepath = "./test_outputs/shape_annealing_test.aj1"
output_filepath = "./test_outputs"
output_filename_no_extension = 'single_objective_output_plots'
new_visualization = VisualizationObject(aj1_filepath=results_filepath, output_filepath=output_filepath,
                                        output_filename_no_extension=output_filename_no_extension)
new_visualization.create_standard_output_single_objective()

results_filepath = "./test_outputs/shape_annealing_multiple_grammar_sets.aj1"
output_filepath = "./test_outputs"
output_filename_no_extension = 'single_objective_multiple_grammar_sets'
new_visualization = VisualizationObject(aj1_filepath=results_filepath, output_filepath=output_filepath,
                                        output_filename_no_extension=output_filename_no_extension)
new_visualization.create_standard_output_single_objective()


# Multi objective output
results_filepath = "./test_outputs/MOSA_test.aj1"
output_filepath = "./test_outputs"
output_filename_no_extension = 'mosa_standard_output_plots'

new_visualization = VisualizationObject(aj1_filepath=results_filepath, output_filepath=output_filepath,
                                        output_filename_no_extension=output_filename_no_extension)
new_visualization.create_standard_output_multi_objective(mesh_or_cylinder='mesh')'''

# Test mass output:
DAED_PATH = "/Users/kodak/Desktop/PERDIX-Mac/DAEDALUS2"
results_filepath = "./test_outputs/MOSA_test.aj1"
output_filepath = "./test_outputs/mass_design_export"
SEQ_filepath = "./M13.txt"

mass_design_export(DAED_PATH, output_filepath, results_filepath, scaffold_sequence_filepath=SEQ_filepath)
print('done')
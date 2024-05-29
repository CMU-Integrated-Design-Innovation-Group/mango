from mango.visualizations.mango_output_visualization import VisualizationObject
from mango.utils.design_io import export_DNA_design
import dill

results_filepath = "/Users/kodak/Desktop/mango/Data Availability/data_output/DNA30/multiobjective_output2.aj1"
output_filepath = "./data_output/plots_files"
output_filename_no_extension = 'multiobjective_output2_DNA30'

new_visualization = VisualizationObject(aj1_filepath=results_filepath, output_filepath=output_filepath,
                                        output_filename_no_extension=output_filename_no_extension)
new_visualization.create_standard_output_multi_objective()

# After creating plots files also export DNA origami designs for oxDNA simulations:
DAED_PATH = '/Users/kodak/Desktop/PERDIX-Mac/DAEDALUS2'
SEQUENCE_FILE = './M13.txt'
with open(results_filepath, 'rb') as f:
    data = dill.load(f)
all_designs = data.MOSA_archive

## NOTE TO SELF: REDO THIS WITH CORRECT DATA POINTS WITH NEW DATA FILE """
'''des1 = all_designs[0.12, 0.815]
des2 = all_designs[0.171, 0.558]
des3 = all_designs[0.343, 0.252]'''
'''
These designs are in the paper:
des1 = all_designs[4.53, 0.1718]
des2 = all_designs[8.43, 0.1548]
des3 = all_designs[15.19, 0.1096]
des4 = all_designs[21.21, 0.0333]
'''
#export_DNA_design(automated_scaffold_executable_path=DAED_PATH, design=des1, export_path='./data_output2',
#                  savename_no_extension='design1_DNA30', scaffold_sequence_filepath=SEQUENCE_FILE)
'''export_DNA_design(automated_scaffold_executable_path=DAED_PATH, design=des2, export_path='./data_output2',
                  savename_no_extension='design2_DNA30', scaffold_sequence_filepath=SEQUENCE_FILE)
export_DNA_design(automated_scaffold_executable_path=DAED_PATH, design=des3, export_path='./data_output2',
                  savename_no_extension='design3_DNA30', scaffold_sequence_filepath=SEQUENCE_FILE)
export_DNA_design(automated_scaffold_executable_path=DAED_PATH, design=des4, export_path='./data_output',
                  savename_no_extension='design4', scaffold_sequence_filepath=SEQUENCE_FILE)'''
print('done')
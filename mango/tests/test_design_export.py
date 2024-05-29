"""
This will show how a user can export specific PLY files for single- and multi-objective problems
"""
import dill
import random
from mango.utils.design_io import export_design_to_ply, export_list_of_designs_to_ply

with open('./test_outputs/shape_annealing_test.aj1', 'rb') as f:
    single_objective_data = dill.load(f)

with open('./test_outputs/MOSA_test.aj1', 'rb') as f:
    multi_objective_data = dill.load(f)


# Export a selected single design (e.g., the final optimized design via shape annealing):
final_design = single_objective_data.design_space
export_design_to_ply(design=final_design,
                     savepath='./test_outputs',
                     savename_no_extension='single_design_output')

# Export a set of selected designs from the archive of a multiobjective problem:
final_archive_designs = list(multi_objective_data.MOSA_archive.keys())
# Selecting 5 random designs to export:
designs_to_export = []
for _ in range(5):
    random_design_key = random.randrange(len(final_archive_designs))
    designs_to_export.append(multi_objective_data.MOSA_archive[final_archive_designs[random_design_key]])

# Call the export:
export_list_of_designs_to_ply(design_list=designs_to_export,
                              savepath='./test_outputs/list_export',
                              list_of_savenames_no_extensions=['Design1', 'Design2', 'Design3', 'Design4', 'Design5'])



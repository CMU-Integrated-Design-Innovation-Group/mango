"""
Automatically create the default plots files for the provided folder paths:
"""
from mango.visualizations.mango_output_visualization import VisualizationObject
import os
from pathlib import Path
from InternalFeatures import edge_cutoff_constraint, porosity_objective, cylinder_volume


fpath = "./data_output/InternalFeatures"
output_filepath = "./data_output/plots_files"


for i in os.listdir(fpath):
    if i != '.DS_Store':
        # Top level folder:
        curPath = os.path.join(fpath, i)
        output_name = Path(i).stem
        new_visualization = VisualizationObject(aj1_filepath=curPath, output_filepath=output_filepath,
                                                output_filename_no_extension=output_name)
        new_visualization.create_standard_output_single_objective(mesh_or_cylinder='cylinder')

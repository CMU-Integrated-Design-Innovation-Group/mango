"""
A.J. Vetturini
IDIG and MMBL
Carnegie Mellon University

This script will support creating an output file that can be used within the kodak visualization tool. This scripts
purpose is to create a simple binary file that can then be read-in by the visualization tool to create the validated
"""
# Import modules
import dill
from dataclasses import dataclass
from mango.optimizers.single_objective_shape_annealing import ShapeAnneal
from mango.visualizations.display_mesh_design import CylindricalRepresentation, MeshRepresentation
import numpy as np
import warnings
import os
from mango.visualizations.optimizer_analysis import MOSA_Analysis

# Import kodak module which is a package I built to visualize my results figures. You can also manually create your
# own results by reading in the results file via dill.load(filepath)
try:
    from kodak.kodak_creator import KodakPlots
except ModuleNotFoundError:
    raise Exception("Missing kodak toolkit package. Please see https://github.com/ajvetturini/kodak_toolkit/tree/master"
                    " for how to install this visualization Python package")


# This returns a list of 33 evenly spaced point (for a total of 35) to use in a design animation for the multiobjective
# output.
def select_evenly_spaced_tuples(tuples_list: list, num_points=33) -> np.array:
    """
    This function finds 33 evenly space points amongst a list of 2D tuples [(1, 2), (3, 4), ...] and uses these points
    as the default exported designs which a user can interact with in the UI window.

    Why 33 as the default? The plots files I create with this value were a reasonable enough size.
    """
    # Always select the first and last tuples
    selected_indices = [0, len(tuples_list) - 1]

    # Calculate the spacing between the tuples to select approximately 'num_points' tuples
    spacing = max(1, len(tuples_list) // (num_points - 1))

    # Select tuples at evenly spaced intervals
    for i in range(spacing, len(tuples_list) - 1, spacing):
        selected_indices.append(i)

    # Extract the selected tuples
    selected_tuples = [tuples_list[i] for i in selected_indices]

    return np.array(selected_tuples)

@dataclass
class VisualizationObject(object):
    """
    This class quickly creates visualization files which can be visualized using the custom kodak viewer found at
    https://ajvetturini.github.io/kodak/. A user can create much more powerful custom scripts, but these defaults
    create a simple, minimal example output that one might be most interested in.

    Parameters
    ------------
    aj1_filepath : String filepath location to the .aj1 results file output from the optimizer
    output_filepath : Where to store the plots file (as only a folderpath)
    output_filename_no_extension : String name of the file to output
    max_number_of_windows : Integer number of windows that the kodak viewer will let you have open. Just here to prevent
                            someone from accidentally opening too many windows, likely don't change the value.
    """
    aj1_filepath: str  # Path to the new AJ1 file to read in
    output_filepath: str
    output_filename_no_extension: str
    max_number_of_windows: int = 10  # The user can change this value, but this is the max # of windows that can be
                                     # created inside of the kodak tool. This may cause memory issues if set too high.

    def __post_init__(self):
        self.data = self.read_in_file()  # Dictionary of data
        if not os.path.exists(self.output_filepath):
            os.mkdir(self.output_filepath)


    def read_in_file(self) -> ShapeAnneal:
        """ Simple function reading in the passed in file and loading the data in via dill """
        with open(self.aj1_filepath, 'rb') as f:
            data = dill.load(f)
        return data


    def create_standard_output_single_objective(self, mesh_or_cylinder: str = 'cylinder') -> None:
        """
        This function will create the standard output windows for a single objective optimization process consisting of:
            - The objective function valuation over the number of epochs which is interactive and a user can select
              from this list
            - A window describing all constraints of the problem

        :param mesh_or_cylinder: Values are either "mesh" or "cylinder" and change the visual appearance of the
                                 generated design
        """
        if mesh_or_cylinder.lower() not in ['mesh', 'cylinder']:
            raise Exception('mesh_or_cylinder can only be "mesh" or "cylinder"')
        kodak_plots = KodakPlots()
        x1, y1 = np.arange(1, len(self.data.objectives_during_annealing), 1), self.data.objectives_during_annealing
        x2, y2 = np.arange(1, len(self.data.random_walk_objective_values), 1), self.data.random_walk_objective_values
        obj1 = kodak_plots.scatter_plot(x=x1, y=y1, lines_markers_both='lines', name='Optimization')
        obj2 = kodak_plots.scatter_plot(x=x2, y=y2, lines_markers_both='lines', name='Random Walk')

        # Add design evolution:
        scatterX, scatterY, scatterXY, stored_designs = [], [], [], []
        for epoch, design_evolution in self.data.design_evolution.items():
            if design_evolution[0] not in y1:
                epoch = 0
            else:
                epoch = y1.index(design_evolution[0])
            scatterX.append(epoch)
            scatterY.append(design_evolution[0])
            scatterXY.append((epoch, design_evolution[0]))
            stored_designs.append(design_evolution[1])

        # First add X Y as scatter points:
        obj3 = kodak_plots.scatter_plot(x=scatterX, y=scatterY, lines_markers_both='markers', name='Design evolution')

        # Next save out the design_evolutions for an interactive plot:
        if len(stored_designs) > 100:
            raise Exception('Too many designs were found in the interactive chart, please reduce this number for'
                            'file size reasons. Create multiple plots files if needed.')
        if len(stored_designs) > 35:
            warnings.warn(
                'This many exported designs will create a large plots file, be aware!')
        print('Exporting designs, depending on length of design this may take a while due to 3D operations')
        design_map = {}

        for xy, des in zip(scatterXY, stored_designs):
            if mesh_or_cylinder == 'cylinder':
                newDes = CylindricalRepresentation(design=des, bounding_box=des.bounding_box)
            else:
                newDes = MeshRepresentation(design=des, bounding_box=des.bounding_box)

            #newDes.create_plot()
            fig = newDes.return_figure()
            design_map[str(xy)] = fig.to_json()
            del newDes

        # Store interactive data:
        kodak_plots.store_interactive_points(design_map=design_map)


        # Store plot data:
        plot_description = 'This plot shows a random walk through space (i.e., by randomly applying grammars) compared ' \
                           'to the optimization process. This helps visualize that the optimization is actually ' \
                           'working. Each of the points is selectable and will "open" a window showing the design at ' \
                           'that stage of the optimization process.'

        kodak_plots.standard_2D_layout.xaxis.title = 'Iteration #'
        kodak_plots.standard_2D_layout.yaxis.title = self.data.objective_function.name
        kodak_plots.add_new_plot(traces=[obj1, obj2, obj3], window_title='Generative process results',
                                 description=plot_description)


        ### Add problem definition:
        if self.data.objective_function.extra_params is None:
            self.data.objective_function.extra_params = 'None'
        kodak_plots.add_problem_definition(objective_function=self.data.objective_function,
                                           design_constraints=self.data.design_constraints)

        ## Write out data:
        kodak_plots.write_json_output(write_directory=self.output_filepath,
                                      savename_no_extension=self.output_filename_no_extension)


    def create_standard_output_multi_objective(self, mesh_or_cylinder: str = 'cylinder') -> None:
        """
        This function will create a standard output for a multiobjective optimization process
            - The found Pareto front which is interactable and allows a user to explore design variation
            - An animation of the pareto growth
            - A window describing all constraints of the problem

        :param mesh_or_cylinder: Values are either "mesh" or "cylinder" and change the visual appearance of the
                                 generated design
        """
        if mesh_or_cylinder.lower() not in ['mesh', 'cylinder']:
            raise Exception('mesh_or_cylinder can only be "mesh" or "cylinder"')
        if hasattr(self.data, 'objective_functions'):
            if len(self.data.objective_functions) > 2:
                raise Exception('Currently this standard output only supports 2-objective output windows. '
                                'You will have to read the documentation about the .aj1 output to create a proper kodak plots output.')
        else:
            raise Exception('Invalid aj1 data file.')
        kodak_plots = KodakPlots()

        if hasattr(self.data, 'MOSA_archive') and hasattr(self.data, 'archive_post_temperature_initialization') and hasattr(self.data, 'objective_functions'):
            final_pareto = np.array(list(self.data.MOSA_archive.keys()))
            random_walk_pareto = np.array(self.data.archive_post_temperature_initialization)
            # Now we will store a random selection of designs evenly spaced along the Pareto. We do this because there
            # are way too many designs in the final pareto that too large a file could be created.
            if len(final_pareto) > 35:
                selected_designs = select_evenly_spaced_tuples(list(final_pareto))
                # In this case, we then want to remove the values of selected_design from final-pareto:
                dtype = [('col1', final_pareto.dtype), ('col2', selected_designs.dtype)]
                struct_array1 = final_pareto.view(dtype).reshape(-1)
                struct_array2 = selected_designs.view(dtype).reshape(-1)
                not_in_array2 = ~np.isin(struct_array1, struct_array2)
                final_pareto = final_pareto[not_in_array2.view(bool).reshape(final_pareto.shape[0])]

            else:
                # If this final pareto was less than the final size we just use whatever the final pareto is
                selected_designs = final_pareto

            selectedX, selectedY = selected_designs[:, 0], selected_designs[:, 1]
            # Next create traces:
            xp, yp = final_pareto[:, 0], final_pareto[:, 1]
            xr, yr = random_walk_pareto[:, 0], random_walk_pareto[:, 1]
            obj1 = kodak_plots.scatter_plot(x=xp, y=yp, lines_markers_both='markers', opacity=0.7, name='Found Pareto',
                                            marker=dict(size=8, color='rgb(68, 170, 153)', symbol='circle',
                                                        line=dict(width=2, color='Black')))
            obj2 = kodak_plots.scatter_plot(x=xr, y=yr, lines_markers_both='markers', name='Random Walk',
                                            marker=dict(size=8, color='rgb(204, 102, 119)', symbol='circle',
                                                        line=dict(width=2, color='Black'))
                                            )
            obj3 = kodak_plots.scatter_plot(x=selectedX, y=selectedY, lines_markers_both='markers',
                                            name='Select Pareto design', marker=dict(size=8, color='rgb(17, 119, 51)',
                                                                                     symbol='diamond',
                                                        line=dict(width=2, color='Black')))

            o1, o2 = self.data.objective_functions[0].name, self.data.objective_functions[1].name
            kodak_plots.standard_2D_layout.xaxis.title = o1
            kodak_plots.standard_2D_layout.yaxis.title = o2


            # Now adding in design evolution as selectable features:
            print('Exporting designs, this may take a few minutes.')
            design_map = {}

            # Loop over the stored designs:
            for des in selected_designs:
                poly_design = self.data.MOSA_archive[tuple(des)]
                if mesh_or_cylinder == 'cylinder':
                    newDes = CylindricalRepresentation(design=poly_design, bounding_box=poly_design.bounding_box)
                else:
                    newDes = MeshRepresentation(design=poly_design, bounding_box=poly_design.bounding_box)
                # newDes.create_plot()
                fig = newDes.return_figure()
                design_map[str(tuple(des))] = fig.to_json()
                del newDes

            # Store interactive data:
            kodak_plots.store_interactive_points(design_map=design_map)

            # Store plot data:
            plot_description = 'This plot shows a sample Pareto found using the mango generative design package. Here' \
                               ' the objective functions are a convexity measure (volume of mesh divided by the volume' \
                               ' of the convex hull of the points in the mesh) and the surface area to volume ratio of' \
                               ' the mesh. The tradeoff is essentially the volume of the mesh which can be explored via' \
                               ' the interactive Pareto points.'
            kodak_plots.add_new_plot(traces=[obj1, obj2, obj3], window_title='Found Pareto front',
                                     description=plot_description)

            animated_figure = self.create_pareto_animation_figure()
            # Add to kodak_plots:
            animation = {
                '_title': 'Found Pareto animation during generative process',
                '_closeable': False,
                '_description': 'Animation showing the stored Pareto points during the generative process',
                '_showGraphSettingsBar': 'scatter3d',
                '_data': animated_figure
            }
            kodak_plots.all_plots['Found Pareto animation during generative process'] = animation


            ### Add problem definition:
            condensed_format = Obj(self.data.objective_functions)
            kodak_plots.add_problem_definition(objective_function=condensed_format,
                                               design_constraints=self.data.design_constraints)

            ## Write out data:
            kodak_plots.write_json_output(write_directory=self.output_filepath,
                                          savename_no_extension=self.output_filename_no_extension)
        else:
            raise Exception('Invalid aj1 data file.')


    def create_pareto_animation_figure(self) -> str:
        """ Simple function used to create the relevant animation figure for the output file """
        mosa_analysis = MOSA_Analysis(results_path=self.aj1_filepath)
        mosa_analysis.create_pareto_animation_figure()
        return mosa_analysis.return_fig().to_json()


@dataclass
class Obj(object):
    """
    This class is not end-user facing, it is used to format an object for input
    """
    objs: list

    def __post_init__(self):
        self.name = ''
        self.extra_params = {}

        for obj in self.objs:
            curName = obj.name
            curParams = obj.extra_params
            if self.name == '':
                self.name = curName
            else:
                self.name = self.name + ' and ' + curName

            if curParams is not None:
                for k, v in curParams.items():
                    self.extra_params[k] = v
        # Now after this if extra_params is still empty we set it to "none"
        if self.extra_params == {}:
            self.extra_params = 'None'
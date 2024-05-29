"""
A.J. Vetturini
IDIG and MMBL
Carnegie Mellon University

This script is responsible for the creation of plots that a user may be interested in creating for their results file.
While many results are created with this script, you may want to create your own plotting functions after the fact once
you know what type of data you want to specifically analyze!

NOTE:
    Some of the code here is HORRIBLE. I intially wrote these internally for a completely different use-case that I
    later scrapped. The practices used here do not make sense for this application, and I really need to update them
    when I find time (will I find the time?).
"""
# Import modules
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import dill
from dataclasses import dataclass
import numpy as np
from itertools import combinations
import pandas as pd

@dataclass
class MOSA_Analysis(object):
    """
    This class creates some simple analysis of the optimizer performance such that a user can readily assess how "well"
    the generative process went.

    Parameters
    ------------
    results_path : Filepath to .aj1 file that conducted a MOSA optimization

    DISPLAY options that control colors / font sizes:
    tick_font_size : Integer point-size of ticker font (default: 16)
    axis_title_font_size : Integer point-size of axis titles (default: 18)
    annotation_font_size : Integer point-size of annotation / text boxes (default: 14)
    sub_plot_title_font_size : Integer point-size of subplot titles (default: 18)
    plot_title_font_size : Integer point-size of figure title (default: 18)
    legend_font : Integer point-size of legend titles (default: 14)
    animation_epoch_step_size : Essentially the result-save frequency (default: 1)
    """
    results_path: str  # This is the path to read in a file
    tick_font_size: int = 16
    axis_title_font_size: int = 18
    annotation_font_size: int = 14
    sub_plot_title_font_size: int = 14
    plot_title_font_size: int = 18
    legend_font: int = 14
    animation_epoch_step_size: int = 1


    def __post_init__(self):
        self.obj_name_function_map = {}
        self.obj_name_ct_map = {}
        self.ct_obj_name_map = {}
        with open(self.results_path, 'rb') as f:
            self.results_data = dill.load(f)

        # After reading in the data, create a dictionary of objective_function_names
        for ct, obj in enumerate(self.results_data.objective_functions):
            self.obj_name_function_map[obj.name] = self.results_data.objective_functions[ct]
            self.obj_name_ct_map[obj.name] = ct
            self.ct_obj_name_map[ct] = obj.name

        self.master_layout = go.Layout(
            autosize=True,
            xaxis=go.layout.XAxis(linecolor='rgba(0, 0, 0, 1)',
                                  linewidth=2,
                                  mirror=True,
                                  showgrid=False,
                                  title_font=dict(size=self.axis_title_font_size, family='Helvetica', color='black'),
                                  tickfont=dict(size=self.tick_font_size, family='Helvetica', color='black'),
                                  ),

            yaxis=go.layout.YAxis(linecolor='rgba(0, 0, 0, 1)',
                                  linewidth=2,
                                  mirror=True,
                                  showgrid=False,
                                  title_font=dict(size=self.axis_title_font_size, family='Helvetica', color='black'),
                                  tickfont=dict(size=self.tick_font_size, family='Helvetica', color='black'),
                                  ),
            plot_bgcolor='white'
        )
        self.master_fig = go.Figure(layout=self.master_layout)


    def objective_function_trace(self, objective_function_name: str) -> go.Scatter:
        """
        This function will use plotly to plot an objective function over the number of iterations the process ran for.

        :param objective_function_name: String name of the objective function being plotted
        :return: A plotly trace of objective function value vs iteration #
        """
        # We will create side-by-side plots so we can see how the pareto optimal surface came to be during MOSA:
        if objective_function_name not in self.obj_name_ct_map:
            raise Exception('Unable to find the name of this function inside of the mapping. Check to see what the '
                            'stored function name is when you ran MOSA.py!')
        idx = self.obj_name_ct_map[objective_function_name]
        x = np.arange(1, len(self.results_data.tracked_objectives_all_temp[idx]) + 1)
        trace = go.Scatter(x=x, y=self.results_data.tracked_objectives_all_temp[idx], mode='lines',
                           line=dict(color="black", dash='solid'),
                           showlegend=True, name=objective_function_name)
        return trace

    @staticmethod
    def create_subplot_grid(length):
        """
        This function determines the number of subplots to use for a given multiobjective optimization result.

        :param length: String name of the objective function being plotted
        :return: Number of rows and columns to use in the subplot
        """
        # Find factors of the length
        if length == 2:
            num_rows = 1
            num_cols = 2
        elif length == 3:
            num_rows = 2
            num_cols = 2
        else:
            # Otherwise, algorithmically define this (WIP):
            factors = [i for i in range(1, length + 1) if length % i == 0]

            # Choose a combination that best fits your requirements
            num_rows = min(factors)
            num_cols = length // num_rows

        return num_rows, num_cols

    def plot_all_objective_func_values(self) -> None:
        """
        This function is called to plot all objective function vs iteration #'s for all N objective functions.
        """
        num_objectives = len(self.results_data.objective_functions)
        num_rows, num_cols = self.create_subplot_grid(num_objectives)
        fig = make_subplots(rows=num_rows, cols=num_cols, print_grid=False)
        for i in range(1, num_objectives + 1):
            row_num = (i - 1) // num_cols + 1
            col_num = (i - 1) % num_cols + 1
            obj_name = self.ct_obj_name_map[i-1]
            trace = self.objective_function_trace(objective_function_name=obj_name)
            fig.add_trace(trace, row=row_num, col=col_num)
            # Update figure axes and things:
            fig.update_xaxes(title_text='Iteration #', row=row_num, col=col_num, linecolor='black',
                             linewidth=1, mirror=True, title_font=dict(size=self.axis_title_font_size),
                             tickfont=dict(size=self.tick_font_size), showgrid=False, gridcolor='white')
            fig.update_yaxes(title_text=obj_name, row=row_num, col=col_num, linecolor='black',
                             linewidth=1, mirror=True, title_font=dict(size=self.axis_title_font_size),
                             tickfont=dict(size=self.tick_font_size), showgrid=False, gridcolor='white')

        fig.update_layout(title='Recorded objective function values during MOSA for all design iterations',
                          title_font=dict(size=self.plot_title_font_size), plot_bgcolor='white')
        fig.update_layout(
            legend_font=dict(size=self.legend_font),
        )
        self.master_fig = fig

    def plot_acceptance_probabilities(self) -> None:
        """
        This function will use plotly to plot (and save) the acceptance probabilty during MOSA to verify that we are
        attempting to only accept the objectively good solutions
        """
        # We will create side-by-side plots so we can see how the pareto optimal surface came to be during MOSA:
        fig = make_subplots(rows=1, cols=2,
                            subplot_titles=("Avg. Acceptance Probability per Epoch", "Actual Probability Iterations"),
                            start_cell="top-left")

        x = np.arange(1, len(self.results_data.acceptance_per_epoch) + 1)
        x2 = np.arange(1, len(self.results_data.p_accept_list) + 1)
        fig.add_trace(go.Scatter(x=x, y=self.results_data.acceptance_per_epoch, mode='lines',
                                 line=dict(color="black", dash='solid'), name='Total Acceptance Probability',
                                 showlegend=True), row=1, col=1)

        fig.add_trace(go.Scatter(x=x2, y=self.results_data.p_accept_list, mode='markers',
                                 marker=dict(color='darkgray', line_width=1, size=10), opacity=0.6,
                                 showlegend=False), row=1, col=2)
        fig.add_trace(go.Scatter(x=x2, y=[1] * len(x2), mode='lines', line=dict(color="red", dash='dash'),
                                 showlegend=False), row=1, col=2)

        fig.update_layout(plot_bgcolor='white')
        # Update axes:
        fig.update_xaxes(title_text='Epoch #', title_font=dict(size=self.axis_title_font_size),
                         tickfont=dict(size=self.tick_font_size), linecolor='black', linewidth=1,
                         mirror=True, row=1, col=1)
        fig.update_yaxes(title_text='"Worse Move" Acceptance probability at Epoch',
                         title_font=dict(size=self.axis_title_font_size), tickfont=dict(size=self.tick_font_size),
                         linecolor='black', linewidth=1, mirror=True, row=1, col=1)

        fig.update_xaxes(title_text='"Worse Move" Iteration #', title_font=dict(size=self.axis_title_font_size),
                         tickfont=dict(size=self.tick_font_size), linecolor='black', linewidth=1, mirror=True,
                         row=1, col=2)
        fig.update_yaxes(title_text='Actual Probability Calculation', title_font=dict(size=self.axis_title_font_size),
                         range=[0, 1.1], tickfont=dict(size=self.tick_font_size), linecolor='black', linewidth=1,
                         mirror=True, row=1, col=2)
        fig.update_annotations(font=dict(size=self.plot_title_font_size))
        self.master_fig = fig

    def temperature_trace(self, objective_function_name: str) -> go.Scatter:
        """
        This function will use plotly to plot the temperature of an objective function vs iteration #

        :param objective_function_name: ObjectiveFunction name attribute (as a string)
        :return: Plotly trace of temperature vs epoch #
        """
        # We will create side-by-side plots so we can see how the pareto optimal surface came to be during MOSA:
        if objective_function_name not in self.obj_name_ct_map:
            raise Exception('Unable to find the name of this function inside of the mapping. Check to see what the '
                            'stored function name is when you ran MOSA.py!')
        idx = self.obj_name_ct_map[objective_function_name]
        x = np.arange(1, len(self.results_data.temp_tracker[idx]) + 1)
        trace = go.Scatter(x=x, y=self.results_data.temp_tracker[idx], mode='lines',
                           line=dict(color="black", dash='solid'),
                           showlegend=True, name=objective_function_name)
        return trace

    def plot_temperatures_over_time(self) -> None:
        """
        This function is called to plot all objective function temperature's vs epoch #'s.
        """
        num_objectives = len(self.results_data.objective_functions)
        num_rows, num_cols = self.create_subplot_grid(num_objectives)
        fig = make_subplots(rows=num_rows, cols=num_cols, print_grid=False,
                            subplot_titles=tuple(self.obj_name_ct_map.keys()))
        for i in range(1, num_objectives + 1):
            row_num = (i - 1) // num_cols + 1
            col_num = (i - 1) % num_cols + 1
            obj_name = self.ct_obj_name_map[i-1]
            trace = self.temperature_trace(objective_function_name=obj_name)
            fig.add_trace(trace, row=row_num, col=col_num)
            # Update figure axes and things:
            fig.update_xaxes(title_text='Epoch #', row=row_num, col=col_num, linecolor='black',
                             linewidth=1, mirror=True, title_font=dict(size=self.axis_title_font_size),
                             tickfont=dict(size=self.tick_font_size), showgrid=False, gridcolor='white')
            fig.update_yaxes(title_text="Temperature", row=row_num, col=col_num, linecolor='black',
                             linewidth=1, mirror=True, title_font=dict(size=self.axis_title_font_size),
                             tickfont=dict(size=self.tick_font_size), showgrid=False, gridcolor='white')

        fig.update_layout(title='Recorded temperature values during MOSA',
                          title_font=dict(size=self.plot_title_font_size), plot_bgcolor='white')
        fig.update_layout(
            legend_font=dict(size=self.legend_font),
        )
        fig.update_annotations(font=dict(size=self.sub_plot_title_font_size))
        self.master_fig = fig


    @staticmethod
    def calculate_pareto_front(df: pd.DataFrame, field1: str, field2: str) -> pd.DataFrame:
        """
        This function calculated the actual pareto front found in the dataframe passed. This is because when 3 or more
        objectives are used, then the points in the archive are not always pareto-optimal points for the current plot!

        :param df: Pandas dataframe containing the obj values of the current archive
        :param field1: First objective function name ("Y" axis of the Pareto)
        :param field2: Second objective function name ("X" axis of the Pareto)
        :return: Pandas dataframe of non-dominating points
        """
        # Initialize a boolean mask to identify non-dominated points
        is_non_dominated = np.ones(len(df), dtype=bool)

        for i, row in df.iterrows():
            if is_non_dominated[i]:
                # Check if there is any other point that dominates the current point
                dominates = (
                        (df[field1].values <= row[field1]) &
                        (df[field2].values <= row[field2]) &
                        ((df[field1].values < row[field1]) | (df[field2].values < row[field2]))
                )
                dominates[i] = False  # Don't compare point to self
                is_non_dominated[dominates] = False  # Mark dominated points

        # Use the mask to filter the dataframe for non-dominated points
        non_dominated_df = df[is_non_dominated]
        return non_dominated_df


    def pareto_2D_section(self, pts: list, field1: str, field2: str, color: str, opac: float, showleg=True) -> go.Scatter:
        """
        This function will create a 2D Pareto trace provided the 2 field names for each function

        :param pts: List of all points that might be in the Pareto
        :param field1: String of first objective function name ("Y" axis of the Pareto)
        :param field2: String of second objective function name ("X" axis of the Pareto)
        :param color: String plotly color or 'rgb(R, G, B)' of the trace to use
        :param opac: Floating point opacity to use for the scatter trace
        :param showleg: Boolean determining if the trace should be shown in the legend (default TRUE)
        :return: Plotly scatter trace of 2D pareto section to plot
        """
        f1_idx, f2_idx = self.obj_name_ct_map[field1], self.obj_name_ct_map[field2]
        new_pts = [tuple(item[index] for index in [f1_idx, f2_idx]) for item in pts]
        df = pd.DataFrame(new_pts, columns=[field1, field2])
        pareto_df = self.calculate_pareto_front(df=df, field1=field1, field2=field2)
        if opac == 1:
            name = 'Final Pareto Front Found'
        else:
            name = 'Initial Pareto Front from Random Walk'
        pareto_df = pareto_df.sort_values(by=field1)  # Sort points for a "clean" line
        pareto_pts = go.Scatter(x=pareto_df[field1], y=pareto_df[field2], mode='markers+lines', name=name,
                                legendgroup='animation', marker=dict(color=color, line_width=1, size=6),
                                opacity=opac, showlegend=showleg)
        return pareto_pts


    def plot_final_2D_paretos(self) -> None:
        """
        This function is used to plot the final pareto plot(s). If there was only 2 specified functions, then 1 figure
        is returned. If 3 functions are used, then 3 plots are created and so on. More functions will return more 2D
        cross sections, be warned! We only recommend MOSA for 2 or 3 objectives in its current implementation.
        """
        num_objectives = len(self.results_data.objective_functions)
        if num_objectives == 2:
            fig = go.Figure()
            f1, f2 = self.results_data.objective_functions[0], self.results_data.objective_functions[1]
            random_walk_pareto = self.results_data.archive_post_temperature_initialization
            final_pareto_idx = list(self.results_data.pareto_animation.keys())
            final_pareto_found = self.results_data.pareto_animation[final_pareto_idx[-1]]
            pareto_initial_points = self.pareto_2D_section(pts=random_walk_pareto, field1=f1.name,
                                                           field2=f2.name, color='red', opac=0.5)
            pareto_final_points = self.pareto_2D_section(pts=final_pareto_found, field1=f1.name,
                                                         field2=f2.name, color='forestgreen', opac=1)
            fig.add_traces(data=[pareto_initial_points, pareto_final_points])
            fig.update_xaxes(title_text=f1.name, linecolor='black',
                             linewidth=1, mirror=True, title_font=dict(size=self.axis_title_font_size),
                             tickfont=dict(size=self.tick_font_size), showgrid=False, gridcolor='white')
            fig.update_yaxes(title_text=f2.name, linecolor='black',
                             linewidth=1, mirror=True, title_font=dict(size=self.axis_title_font_size),
                             tickfont=dict(size=self.tick_font_size), showgrid=False, gridcolor='white')


        elif num_objectives == 3:
            num_rows, num_cols = self.create_subplot_grid(num_objectives)
            fig = make_subplots(rows=num_rows, cols=num_cols, print_grid=False)
            i = 1
            showleg = True
            for f1, f2 in combinations(self.results_data.objective_functions, 2):
                row_num = (i - 1) // num_cols + 1
                col_num = (i - 1) % num_cols + 1
                random_walk_pareto = self.results_data.archive_post_temperature_initialization
                final_pareto_idx = list(self.results_data.pareto_animation.keys())
                final_pareto_found = self.results_data.pareto_animation[final_pareto_idx[-1]]
                pareto_initial_points = self.pareto_2D_section(pts=random_walk_pareto, field1=f1.name, showleg=showleg,
                                                               field2=f2.name, color='red', opac=0.5)
                pareto_final_points = self.pareto_2D_section(pts=final_pareto_found, field1=f1.name, showleg=showleg,
                                                             field2=f2.name, color='forestgreen', opac=1)
                # Now we add these ponints to the proper row n column:
                fig.add_trace(pareto_initial_points, row=row_num, col=col_num)
                fig.add_trace(pareto_final_points, row=row_num, col=col_num)
                fig.update_xaxes(title_text=f1.name, linecolor='black', row=row_num, col=col_num,
                                 linewidth=1, mirror=True, title_font=dict(size=self.axis_title_font_size),
                                 tickfont=dict(size=self.tick_font_size), showgrid=False, gridcolor='white')
                fig.update_yaxes(title_text=f2.name, linecolor='black', row=row_num, col=col_num,
                                 linewidth=1, mirror=True, title_font=dict(size=self.axis_title_font_size),
                                 tickfont=dict(size=self.tick_font_size), showgrid=False, gridcolor='white')
                showleg = False
                i += 1

        else:
            raise Exception('We only recommend MOSA for 2 or 3 objectives. If you want plots for higher dimensionality,'
                            ' we recommend you create your own plotting figures as there will be too many 2D plots to'
                            ' analyze the results efficiently!')

        # Depending on "return_or_show" we may return the figure:
        fig.update_layout(title='Pareto optimal front found compared to the random walk pareto',
                          title_font=dict(size=self.plot_title_font_size), plot_bgcolor='white',
                          legend_font=dict(size=self.legend_font))
        self.master_fig = fig


    def plot_final_3D_pareto(self) -> None:
        """
        This function is only used when 3 objectives are used. It will show the surface plot of the random walk pareto
        compared to the final pareto found
        """
        if len(self.results_data.objective_functions) != 3:
            raise Exception('This function is only valid for 3 objective functions in the pareto. ')
        final_archive = list(self.results_data.MOSA_archive.keys())
        initial_archive = self.results_data.archive_post_temperature_initialization
        xf, yf, zf = zip(*final_archive)
        xi, yi, zi = zip(*initial_archive)
        ## SORTING:
        sorted_final = np.argsort(xf)
        sorted_initial = np.argsort(xi)
        xf, yf, zf = np.array(xf)[sorted_final], np.array(yf)[sorted_final], np.array(zf)[sorted_final]
        xi, yi, zi = np.array(xi)[sorted_initial], np.array(yi)[sorted_initial], np.array(zi)[sorted_initial]
        t1 = go.Scatter3d(x=xi, y=yi, z=zi, mode='markers+lines', marker=dict(size=8, color='red', symbol='square'),
                          opacity=0.4, showlegend=True, legendgroup='1', name='Initial Pareto Front from Random Walk')
        t2 = go.Scatter3d(x=xf, y=yf, z=zf, mode='markers+lines', marker=dict(size=8, color='forestgreen', symbol='square'),
                          opacity=1.0, showlegend=True, legendgroup='2', name='Final Pareto Front Found')

        # Update layout and things:
        f1, f2, f3 = tuple(self.results_data.objective_functions)
        layout = go.Layout(
            scene=dict(
                xaxis=dict(
                    title=f1.name,
                    titlefont=dict(size=self.axis_title_font_size, color='black'),  # Customize title font
                    tickfont=dict(size=self.tick_font_size-4, color='black'),  # Customize tick font
                ),
                yaxis=dict(
                    title=f2.name,
                    titlefont=dict(size=self.axis_title_font_size, color='black'),  # Customize title font
                    tickfont=dict(size=self.tick_font_size-4, color='black'),  # Customize tick font
                ),
                zaxis=dict(
                    title=f3.name,
                    titlefont=dict(size=self.axis_title_font_size, color='black'),  # Customize title font
                    tickfont=dict(size=self.tick_font_size-4, color='black'),  # Customize tick font
                ),
            ),
            plot_bgcolor='white'
        )
        # Create a 3D surface plot
        fig = go.Figure(data=[t1, t2], layout=layout)
        fig.update_layout(title='Final 3D Pareto Front Found',
                          title_font=dict(size=self.plot_title_font_size), plot_bgcolor='white',
                          legend_font=dict(size=self.legend_font))
        self.master_fig = fig


    def create_pareto_data(self, pts: list, num_rows: int, num_cols: int, epoch: int, last: bool = False) -> list:
        # Create "initial" frame by population individual plots:
        plot_map = {(0, 0): (0, 1), (0, 1): (1, 2), (1, 0): (0, 2)}
        # Starting data in figure:
        return_trace = []
        for row in range(num_rows):
            for col in range(num_cols):
                if row == (num_rows - 1) and col == (num_cols - 1):
                    # If we are in the final cell, we update the Annotation block
                    if epoch == 0:
                        text_string = f'Epoch: {epoch} (START)'
                    elif not last:
                        text_string = f'Epoch: {epoch}'
                    else:
                        text_string = f'Epoch: {epoch} (END)'
                    return_trace.append((go.Scatter(x=[0.5], y=[0.5], text=[text_string],
                                                    mode="text", showlegend=False,
                                                    textfont=dict(size=self.plot_title_font_size, family='Helvetica', color='black')), (row + 1, col + 1)))
                else:
                    # Otherwise we plot
                    func_to_plot = plot_map[(row, col)]
                    f1, f2 = self.ct_obj_name_map[func_to_plot[0]], self.ct_obj_name_map[func_to_plot[1]]
                    scatter_trace = self.pareto_2D_section(pts=pts, field1=f1, field2=f2, color='red',
                                                           opac=1, showleg=False)
                    # Add to the plot:
                    return_trace.append((scatter_trace, (row + 1, col + 1)))
        return return_trace


    def create_pareto_animation_figure(self) -> None:
        """
        This function initializes + creates the pareto animation figures to see how the Pareto grew over time
        """
        num_objectives = len(self.results_data.objective_functions)
        if num_objectives == 2:
            num_rows, num_cols = 1, 2
        elif num_objectives == 3:
            num_rows, num_cols = 2, 2
        else:
            raise Exception('The animation figure can only be created using 2 or 3 objective functions')

        fig = make_subplots(rows=num_rows, cols=num_cols, print_grid=False)
        fig.update_xaxes(visible=False, row=num_rows, col=num_cols)  # hide the axis and things
        fig.update_yaxes(visible=False, row=num_rows, col=num_cols)
        fig.update_layout(plot_bgcolor='white')
        plot_map = {(0, 0): (0, 1), (0, 1): (1, 2), (1, 0): (0, 2)}
        # Initialize figure layout options:
        for row in range(num_rows):
            if row == num_rows - 1:
                mc = num_cols - 1
            else:
                mc = num_cols
            for col in range(mc):
                func_to_plot = plot_map[(row, col)]
                f1, f2 = self.ct_obj_name_map[func_to_plot[0]], self.ct_obj_name_map[func_to_plot[1]]
                max_f1 = max(self.results_data.tracked_objectives_all_temp[func_to_plot[0]])
                min_f1 = min(self.results_data.tracked_objectives_all_temp[func_to_plot[0]])
                max_f2 = max(self.results_data.tracked_objectives_all_temp[func_to_plot[1]])
                min_f2 = min(self.results_data.tracked_objectives_all_temp[func_to_plot[1]])
                padding = (max_f1 - min_f1) * (5 / 100.0)

                fig.update_xaxes(title_text=f'{f1}', row=row + 1, col=col + 1, linecolor='black', linewidth=2,
                                 mirror=True, title_font=dict(size=self.axis_title_font_size, family='Helvetica', color='black'),
                                 tickfont=dict(size=self.tick_font_size, family='Helvetica', color='black'), range=[min_f1 - padding,
                                                                                 max_f1 + padding])

                fig.update_yaxes(title_text=f'{f2}', row=row + 1, col=col + 1, linecolor='black', linewidth=2,
                                 mirror=True, title_font=dict(size=self.axis_title_font_size, family='Helvetica', color='black'),
                                 tickfont=dict(size=self.tick_font_size, family='Helvetica', color='black'),
                                 range=[min_f2 - padding, max_f2 + padding])
        random_walk_pareto = self.results_data.archive_post_temperature_initialization
        initial_data = self.create_pareto_data(pts=random_walk_pareto, num_rows=num_rows, num_cols=num_cols, epoch=0)
        for data, row_col in initial_data:
            fig.add_trace(data, row=row_col[0], col=row_col[1])


        # Loop over all epochs using the step variable:
        all_frames = []
        last = False
        for epoch in range(1, len(self.results_data.pareto_animation.keys()) + 1, self.animation_epoch_step_size):
            # With each epoch, we create a new dataframe:
            frame = {"data": [], "name": str(epoch)}
            data_to_plot = self.results_data.pareto_animation[epoch]
            if epoch == len(self.results_data.pareto_animation.keys()):
                last = True
            n = self.create_pareto_data(pts=data_to_plot, num_rows=num_rows, num_cols=num_cols, epoch=epoch, last=last)
            data_list = []
            # Add data to frame:
            for d, _ in n:
                data_list.append(d)
            frame["data"] = data_list
            all_frames.append(frame)

        # Ensure last element is processed to the animation:
        if (len(self.results_data.pareto_animation.keys()) - 1) % self.animation_epoch_step_size != 0:
            epoch = list(self.results_data.pareto_animation.keys())[-1] # Grab last index
            frame = {"data": [], "name": str(epoch)}
            data_to_plot = self.results_data.pareto_animation[epoch]
            new_data = self.create_pareto_data(pts=data_to_plot, num_rows=num_rows, num_cols=num_cols, epoch=epoch,
                                               last=True)
            data_list = []
            # Add data to frame:
            for d, _ in new_data:
                data_list.append(d)
            frame["data"] = data_list
            all_frames.append(frame)

        # Finally add frames and slider:
        fig["frames"] = all_frames
        fig["layout"]["updatemenus"] = [
            {
                "buttons": [
                    {
                        "args": [None, {"frame": {"duration": 500, "redraw": False},
                                        "fromcurrent": True, "transition": {"duration": 300,
                                                                            "easing": "quadratic-in-out"}}],
                        "label": "Play",
                        "method": "animate"
                    },
                    {
                        "args": [[None], {"frame": {"duration": 0, "redraw": False},
                                          "mode": "immediate",
                                          "transition": {"duration": 0}}],
                        "label": "Pause",
                        "method": "animate"
                    }
                ],
                "direction": "left",
                "pad": {"r": 10, "t": 87},
                "showactive": False,
                "type": "buttons",
                "x": 0.1,
                "xanchor": "right",
                "y": 0,
                "yanchor": "top"
            }
        ]

        # Depending on "return_or_show" we may return the figure:
        self.master_fig = fig

    def show_fig(self) -> None:
        self.master_fig.show()

    def return_fig(self) -> go.Figure:
        return self.master_fig


@dataclass
class ShapeAnneal_Analysis(object):
    """
    This class creates some simple analysis of the optimizer performance such that a user can readily assess how "well"
    the generative process went.

    Parameters
    ------------
    results_path : Filepath to .aj1 file that conducted a shape annealing optimization

    DISPLAY options that control colors / font sizes:
    tick_font_size : Integer point-size of ticker font (default: 16)
    axis_title_font_size : Integer point-size of axis titles (default: 18)
    annotation_font_size : Integer point-size of annotation / text boxes (default: 14)
    sub_plot_title_font_size : Integer point-size of subplot titles (default: 18)
    plot_title_font_size : Integer point-size of figure title (default: 18)
    legend_font : Integer point-size of legend titles (default: 14)
    """
    results_path: str  # This is the path to read in a file
    tick_font_size: int = 16
    axis_title_font_size: int = 18
    annotation_font_size: int = 14
    sub_plot_title_font_size: int = 14
    plot_title_font_size: int = 18
    legend_font: int = 14


    def __post_init__(self):
        with open(self.results_path, 'rb') as f:
            self.results_data = dill.load(f)

        self.master_layout = go.Layout(
            autosize=True,
            xaxis=go.layout.XAxis(linecolor='rgba(0, 0, 0, 1)',
                                  linewidth=2,
                                  mirror=True,
                                  showgrid=False,
                                  title_font=dict(size=self.axis_title_font_size, family='Arial', color='black'),
                                  tickfont=dict(size=self.tick_font_size, family='Arial', color='black'),
                                  ),

            yaxis=go.layout.YAxis(linecolor='rgba(0, 0, 0, 1)',
                                  linewidth=2,
                                  mirror=True,
                                  showgrid=False,
                                  title_font=dict(size=self.axis_title_font_size, family='Arial', color='black'),
                                  tickfont=dict(size=self.tick_font_size, family='Arial', color='black'),
                                  ),
            plot_bgcolor='white'
        )
        self.master_fig = go.Figure(layout=self.master_layout)

    def plot_objective_trace(self) -> None:
        """
        This function will use plotly to plot (and save) the objective function value during the optimization process.
        """
        # We will create side-by-side plots so we can see how the pareto optimal surface came to be during MOSA:
        x = np.arange(1, len(self.results_data.objectives_during_annealing) + 1)
        trace = go.Scatter(x=x, y=self.results_data.objectives_during_annealing, mode='lines',
                           line=dict(color="black", dash='solid'),
                           showlegend=True, name=self.results_data.objective_function.name)
        self.master_fig = go.Figure(data=[trace], layout=self.master_layout)
        self.master_fig.update_xaxes(title_text='Iteration #')
        self.master_fig.update_yaxes(title_text='Objective Function Value')

    def plot_temperature_trace(self) -> None:
        """
        This function will use plotly to plot (and save) the temperature during the optimization process.
        """
        x = np.arange(1, len(self.results_data.temperature_values) + 1)
        trace = go.Scatter(x=x, y=self.results_data.temperature_values, mode='lines',
                           line=dict(color="black", dash='solid'),
                           showlegend=True, name=self.results_data.objective_function.name)
        self.master_fig = go.Figure(data=[trace], layout=self.master_layout)
        self.master_fig.update_xaxes(title_text='Epoch #')
        self.master_fig.update_yaxes(title_text='Temperature')

    def show_fig(self) -> None:
        self.master_fig.show()

    def return_fig(self) -> go.Figure:
        return self.master_fig


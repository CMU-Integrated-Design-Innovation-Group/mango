"""
This script will allow users to create / edit JSON files for the kodak visualization output tool.

Note: Not all plot types are yet supported and this is actively being developed (thus why this project is not on PyPi
      yet).
"""
import os
import json
import warnings
from typing import Optional, Union
import numpy as np
from dataclasses import dataclass
from kodak.trace_types import PlotTrace, scatter, scatter3d
import plotly.graph_objs as go


@dataclass
class KodakPlots(object):
    # Default parameters a user can change:
    font_style = 'Helvetica'  # Options include: Open Sans, veerdana, arial, sans-serif (probably others, check plotly)
    default_font_color = 'black'

    """
    These methods create the "default" layout for the various plots that Kodak supports using a plotly layout.
    Note: These methods are the "backend" of Kodak because they enable an Academic-style default layout using Plotly
    """
    def create_2d_layout(self, linecolor: str = 'rgba(0, 0, 0, 1)', tick_font_size: int = 16, title_font_size: int = 18,
                         linewidth: float = 2., mirror_line: bool = True, showgrid: bool = False, gridcolor: str = 'black',
                         zeroline: bool = False, plot_bgcolor='rgba(0, 0, 0, 0)', paper_bgcolor='rgba(255,255,255, 1)'):

        return go.Layout(
            height=318,
            width=521,
            xaxis=dict(
                linecolor=linecolor,
                linewidth=linewidth,
                mirror=mirror_line,
                tickfont=dict(size=tick_font_size, family=self.font_style, color=self.default_font_color),
                titlefont=dict(size=title_font_size, family=self.font_style, color=self.default_font_color),
                showgrid=showgrid,
                gridcolor=gridcolor,
                zeroline=zeroline,
            ),
            yaxis=dict(
                linecolor=linecolor,
                linewidth=linewidth,
                mirror=mirror_line,
                tickfont=dict(size=tick_font_size, family=self.font_style, color=self.default_font_color),
                titlefont=dict(size=title_font_size, family=self.font_style, color=self.default_font_color),
                showgrid=showgrid,
                gridcolor=gridcolor,
                zeroline=zeroline,
            ),
            font=dict(family=self.font_style, color='black', size=tick_font_size),
            plot_bgcolor=plot_bgcolor,
            paper_bgcolor=paper_bgcolor,
            autosize=False
        )


    def create_3d_layout(self, display_background: bool = False, tick_font_size: int = 16, axis_font_size: int = 18,
                         backgroundcolor: str = 'rgba(0, 0, 0, 0)', plot_bgcolor: str = 'rgba(0, 0, 0, 0)',
                         paper_bgcolor: str = 'rgba(0, 0, 0, 0)', showticklabels: bool = False, showgrid: bool = False,
                         gridcolor: str = 'black', zerolinewidth: float = 0.):
        return go.Layout(
            scene=dict(
                xaxis_visible=display_background,
                yaxis_visible=display_background,
                zaxis_visible=display_background,
                xaxis=dict(
                    title="X",
                    tickfont=dict(size=tick_font_size, family=self.font_style, color=self.default_font_color),
                    titlefont=dict(size=axis_font_size, family=self.font_style, color=self.default_font_color),
                    showbackground=display_background,
                    showticklabels=showticklabels,
                    showgrid=showgrid,
                    gridcolor=gridcolor,
                    backgroundcolor=backgroundcolor,
                    zerolinewidth=zerolinewidth,

                ),
                yaxis=dict(
                    title="Y",
                    tickfont=dict(size=tick_font_size, family=self.font_style, color=self.default_font_color),
                    titlefont=dict(size=axis_font_size, family=self.font_style, color=self.default_font_color),
                    showbackground=display_background,
                    showgrid=showgrid,
                    showticklabels=showticklabels,
                    gridcolor=gridcolor,
                    backgroundcolor=backgroundcolor,
                    zerolinewidth=zerolinewidth,
                ),
                zaxis=dict(
                    title="Z",
                    tickfont=dict(size=tick_font_size, family=self.font_style, color=self.default_font_color),
                    titlefont=dict(size=axis_font_size, family=self.font_style, color=self.default_font_color),
                    showbackground=display_background,
                    showticklabels=showticklabels,
                    showgrid=showgrid,
                    gridcolor=gridcolor,
                    backgroundcolor=backgroundcolor,
                    zerolinewidth=zerolinewidth,
                ),
            ),
            plot_bgcolor=plot_bgcolor,
            paper_bgcolor=paper_bgcolor,
            margin=dict(l=0, r=0, b=0, t=0)
        )


    def __post_init__(self) -> None:
        self.all_plots = {}  # Dictionary that gets written out
        self.id_count = -1  # Initialize interactive design chart count to -1 since get_id_count will bring to 0

        self.currently_supported = ['scatter', 'scatter3d']  # Here is a running list of supported plotly plot types in
                                                             # kodak. These are the plot.type value.

        ## Creating standard default layouts for all plot styles:
        self.standard_2D_layout = self.create_2d_layout()
        self.standard_3D_layout = self.create_3d_layout()


    def get_id_count(self) -> str:
        self.id_count += 1
        return str(self.id_count)


    """
    This section of code below is for creating plotly figures through my own "wrapper" of sorts to enable the "default"
    academic style feel to the plots. 
    """

    def get_color_and_symbol(self, color_idx: int, symbol_idx: int, symbol_or_dash: str) -> tuple:
        """
        Based on the trace #, a different color and symbol may be used. The user can more clearly "fix" this in the UI
        if a different combo is desired.
        """
        if color_idx > 8:
            warnings.warn('Note: You have many traces in one figure so re-using colors. Reminder that this is an '
                          'academic-minded package, do you want more than 8 traces on one plot?')
            color_idx = color_idx % 8  # We repeat colors
        tol_colors = {  # Color-blind safe color map from: https://personal.sron.nl/~pault/
            1: 'rgb(136, 204, 238)',  # Cyan
            2: 'rgb(68, 170, 153)',  # Teal
            3: 'rgb(17, 119, 51)',  # Green
            4: 'rgb(153, 153, 51)',  # Olive
            5: 'rgb(221, 204, 119)',  # Sand
            6: 'rgb(204, 102, 119)',  # Rose
            7: 'rgb(136, 34, 85)',  # Wine
            8: 'rgb(170, 68, 153)',  # Purple
        }
        if symbol_or_dash == 'symbol':
            if symbol_idx > 5:
                symbol_idx = symbol_idx % 5  # We repeat colors
            symbol_map = {
                1: 'circle',
                2: 'square',
                3: 'diamond',
                4: 'x',
                5: 'cross'
            }
        else:
            if symbol_idx > 4:
                symbol_idx = symbol_idx % 4 # We repeat colors
            symbol_map = {
                1: 'solid',
                2: 'dash',
                3: 'dot',
                4: 'dashdot',
            }

        return tol_colors[color_idx], symbol_map[symbol_idx]


    def scatter_plot(self, x: Optional[Union[list, np.ndarray]], y: Optional[Union[list, np.ndarray]],
                     lines_markers_both: str, **kwargs) -> PlotTrace:
        # First we create a scatter trace based on X and Y:
        if lines_markers_both not in ['lines', 'markers', 'both']:
            raise Exception("Invalid value for lines_markers_both, only options are: 'lines', 'markers', or 'both'")
        marker = go.scatter.Marker()
        line = go.scatter.Line()

        if 'marker' in kwargs.keys():
            marker = kwargs['marker']
            del kwargs['marker']  # We also get rid of it since we are going to pass it in manually
        elif 'marker' not in kwargs.keys() and (lines_markers_both == 'markers' or lines_markers_both == 'both'):
            # NOTE: These values for marker and line get updated later during the figure creation, we need to know
            #       all traces in a figure before determining colors and symbols:
            pass
            #marker = dict(size=8, color='blue', symbol='circle', line=dict(width=2, color='DarkSlateGrey'))

        if 'line' in kwargs.keys():
            marker = kwargs['line']
            del kwargs['line']  # We also get rid of it since we are going to pass it in manually
        elif 'line' not in kwargs.keys() and (lines_markers_both == 'lines' or lines_markers_both == 'both'):
            #line = dict(width=2, color='blue', dash='solid')
            pass

        if lines_markers_both == 'lines':
            return scatter(x=x, y=y, lines_markers_both=lines_markers_both, line=line, **kwargs)
        elif lines_markers_both == 'markers':
            return scatter(x=x, y=y, lines_markers_both=lines_markers_both, marker=marker, **kwargs)
        else:
            return scatter(x=x, y=y, lines_markers_both=lines_markers_both, marker=marker, line=line, **kwargs)


    def scatter3d_plot(self, x: Optional[Union[list, np.ndarray]], y: Optional[Union[list, np.ndarray]],
                       **kwargs) -> PlotTrace:
        return scatter3d(x=x, y=y, **kwargs)


    """
    This section of code is responsible for creating the actual individual plot_data (or WindowObjects) 
    """

    def validate_traces(self, traces: list) -> tuple:
        traces_2D, traces_3D = False, False
        for trace in traces:
            if not isinstance(trace, PlotTrace):
                return False, f'Invalid plotly trace provided which has no "type" attribute: {trace}'

            trace_type = trace.type
            if '3d' in trace_type:
                if traces_2D:
                    return False, 'Can not combine 2D and 3D traces on the same plot.'
                traces_3D = True

            else:
                if traces_3D:
                    return False, 'Can not combine 2D and 3D traces on the same plot.'
                traces_2D = True

            # If we did not error-out, we finally just check if this trace_type is supported:
            if trace_type not in self.currently_supported:
                return False, f'Since this repo is under active development, not all plots are currently supported.' \
                              f'Currently, the {trace_type} plot type is not yet supported. Please open a request!'

        # If we make it here, we return true:
        return True, ''


    def read_traces(self, traces: list) -> tuple:
        """
        This function will read in the list of PlotTraces that the user created and create the proper list and extract
        the required information. As we add more traces, need to add further logic to this function.
        :param traces: List of PlotTraces from trace_types
        :return: List of extracted plotly traces for this figure, the default layout to use, and the plot_type
        """
        extracted_traces = []
        plot_type = None
        color_count, symbol_count = 1, 1
        for trace in traces:
            # First update plot_type:
            if plot_type is None:
                plot_type = trace.type
            elif trace.type not in plot_type:
                plot_type += f'+{trace.type}'

            plt = trace.plotly_trace  # Extract this as we need to update the color n symbol n things:
            cur_mode = plt.mode
            # If either of these is not None, then we know that we want to update their marker or line style:
            if plt.marker.color is None and plt.line.color is None:
                if cur_mode in ['markers', 'markers+lines', 'lines+markers']:
                    color, symb = self.get_color_and_symbol(color_idx=color_count, symbol_idx=symbol_count,
                                                            symbol_or_dash='symbol')
                    new_marker = dict(size=8, color=color, symbol=symb, line=dict(width=2, color='DarkSlateGrey'))
                    trace.plotly_trace.marker = new_marker
                elif cur_mode in ['lines', 'markers+lines', 'lines+markers']:
                    color, symb = self.get_color_and_symbol(color_idx=color_count, symbol_idx=symbol_count,
                                                            symbol_or_dash='dash')
                    new_line = dict(width=2, color=color, dash=symb)
                    trace.plotly_trace.line = new_line

            elif plt.marker.color is None:
                color, symb = self.get_color_and_symbol(color_idx=color_count, symbol_idx=symbol_count,
                                                        symbol_or_dash='symbol')
                new_marker = dict(size=8, color=color, symbol=symb, line=dict(width=2, color='DarkSlateGrey'))
                trace.plotly_trace.marker = new_marker

            elif plt.line.color is None:
                color, symb = self.get_color_and_symbol(color_idx=color_count, symbol_idx=symbol_count,
                                                        symbol_or_dash='dash')
                new_line = dict(width=2, color=color, dash=symb)
                trace.plotly_trace.line = new_line

            # Next extract the traces and update the counts
            extracted_traces.append(trace.plotly_trace)
            color_count += 1
            symbol_count += 1

        # Now determine the default layout. Add logic here for when adding compatibility for "mix and match" plots such
        # as scatter with bar chart:
        if '3d' in plot_type:
            layout = self.standard_3D_layout
        else:
            layout = self.standard_2D_layout

        return extracted_traces, layout, plot_type


    def add_new_plot(self, traces: list, window_title: str, description: str = None, closeable: bool = False,
                           custom_layout: go.Layout = None) -> None:
        """
        This function specifically iterates over the plots_list and calls the proper functions to "add" this data
        to the all_plots dictionary.

        :param traces: All traces located on the plot
        :param window_title: What you want the window title to be
        :param description: The description of the plot
        :param closeable: If the plot should be "closeable"
        :param custom_layout: If the user has a custom layout they may want to use, then pass it in here.
        :return: Nothing, this will just add the plot to the all_plots dictionary:
        """
        # First validate traces can be added to one plot:
        valid, message = self.validate_traces(traces=traces)

        if valid:
            extracted_traces, layout, plot_type = self.read_traces(traces=traces)
            # If user specified a specific layout to use, then we just use that instead:
            if custom_layout:
                layout = custom_layout

            # Finally, we add the new plot to all_plots for export:
            fig = go.Figure(data=extracted_traces, layout=layout)
            self.all_plots[window_title] = {
                '_title': window_title,
                '_closeable': closeable,
                '_description': description if description is not None else 'No plot description provided.',
                '_showGraphSettingsBar': plot_type,
                '_data': fig.to_json()  # Convert to JSON so that I can read it in
            }
        else:
            raise Exception(f'Invalid list of traces passed into function. Error: {message}')


    """
    Below are some custom-made functions that I use Kodak for when creating visualizations from my mango project. 
    
    Since these visualizations are technically broadly applicable, i added the code to this project. If you have any 
    questions / want to more widely apply the code below to different functionalities, please email me or open a PR.
    (avetturi@andrew.cmu.edu)
    """

    def add_problem_definition(self, objective_function, design_constraints):
        problem_definition = {
            'objective_name': objective_function.name,
            'objective_constants': objective_function.extra_params,
            'constraint_edge_lengths': [design_constraints.min_edge_length, design_constraints.max_edge_length],
            'constraint_min_angle': design_constraints.min_face_angle,
            'constraint_max_basepairs': design_constraints.max_number_basepairs_in_scaffold,
        }
        if design_constraints.extra_constraints is not None:
            for ct, extra_constraint in enumerate(design_constraints.extra_constraints):
                f1 = 'custom_constraint_' + str(ct) + '_name'
                f2 = 'custom_constraint_' + str(ct) + '_params'
                problem_definition[f1] = extra_constraint.name
                problem_definition[f2] = extra_constraint.extra_params

        self.all_plots['problem_definition'] = {
            '_title': "Problem specified to mango framework",
            '_closeable': False,
            '_description': 'This figure shows an the various parameters used in the optimization process. A more '
                            'comprehensive list can be seen using the advanced output creation.',
            '_showGraphSettingsBar': 'None',
            '_data': problem_definition
        }


    def store_interactive_points(self, design_map: dict) -> None:
        # This function simply stores mapped data to an interactive_design:
        self.all_plots[f'INTERACTIVE_DESIGN_{self.get_id_count()}'] = design_map

    """
    I/O Functionalities
    """
    def write_json_output(self, write_directory: str, savename_no_extension: str) -> None:
        sname = savename_no_extension + '.plots'
        full_save_path = os.path.join(write_directory, sname)
        with open(full_save_path, 'w') as f:
            json.dump(self.all_plots, f, indent=4)


    def read_json_input(self, filepath: str) -> None:
        # Note: Filepath should be the full path to the file we are reading-in
        with open(filepath, 'r') as f:
            self.all_plots = json.load(f)

        # After reading-in, we need to almost do the "inverse" of __post_init__ to set values properly:
        fin = -1  # Default value
        for plot in self.all_plots.keys():
            if 'INTERACTIVE_DESIGN' in plot:
                fin = plot.split('_')[2]

        self.id_count = fin


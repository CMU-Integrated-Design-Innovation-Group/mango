"""
A.J. Vetturini
IDIG and MMBL
Carnegie Mellon University

This is a visualization function that will utilize plotly to show a design as a mesh file both as a cylindrical
representation and as a standard mesh file
"""
# Import Modules
from mango.design_spaces.bounding_box import TriclinicBox, MonoclinicBox, OrthorhombicBox, TetragonalBox, RhombohedralBox, HexagonalBox, CubicBox
from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
from dataclasses import dataclass
from mango.utils.mango_math import length_and_direction_between_nodes, rotation_matrix_from_vectors, find_smaller_z_value
import plotly.graph_objs as go
import numpy as np
from typing import Union

@dataclass
class MeshRepresentation(object):
    """
    This class lets a user visualize the resultant structure using plotly traces. This function is also called by the
    plots / output file creator, but it has many options that can be adjusted as discussed in the below parameters. This
    visualization is built using plotly and the figures will all show in a web browser.

    Parameters
    ------------
    design : Design we are visualizing of the PolyhedralSpace data object (see design_spaces.polyhedral_design_space).
    bounding_box : The bounding box defined, can be accessed by the design.bounding_box

    BOOLEAN options that turn on / off various visualizations:
    display_node_names : Determines if node names are shown, used for debugging (default: FALSE)
    show_bounding_box : Determines if bounding box is shown, used for debugging (default: TRUE)
    show_mesh_nodes : Determines if nodes are shown as scatter points or not (default: TRUE)
    show_excluded_regions : Determines if the excluded regions are shown or not (default: TRUE)
    opaque_excluded_region : Determines if excluded regions are opaque or not (default: FALSE / not opaque)
    display_background : Determines if a background color is to be used (default: False / blank background)
    display_axes : Determines if the axes and ticker labels are shown (default: True)


    DISPLAY options that control colors / font sizes:
    face_color : String color of mesh faces defined by either plotly color or 'rgb(R, G, B)' (default: darkgray)
    constant_edge_color : String color of edges shown by either plotly color or 'rgb(R, G, B)' (default: None / black)
    symbol_size : Integer size of mesh nodes (default: 6)
    line_width : Integer width of edges shown (default: 4)
    tick_font : Integer point-size of ticker font (default: 14)
    axis_title_font : Integer point-size of axis titles
    camera_position : Camera dictionary (see plotly) for holding the camera angle constant
    node_color :  String color of nodes defined by either plotly color or 'rgb(R, G, B)' (default: green)
    symbol : String symbol of nodes defined by plotly (default: square)
    excluded_color : String color of excluded regions by either plotly color or 'rgb(R, G, B)' (default: red)
    plot_background_color : String color of background defined by either plotly color or 'rgb(R, G, B)' (default: None)
    """
    design: PolyhedralSpace
    bounding_box: Union[TriclinicBox, MonoclinicBox, OrthorhombicBox, TetragonalBox,
                        RhombohedralBox, HexagonalBox, CubicBox]
    display_node_names: bool = False
    show_bounding_box: bool = True
    show_mesh_nodes: bool = True
    show_excluded_regions: bool = True
    display_background: bool = False  # Determine if axis titles and ticker labels are shown
    display_axes: bool = True  # Determine if axis titles and ticker labels are shown
    face_color: str = None  # Override potential face color of mesh if desired
    constant_edge_color: str = None  # If set then all edges will appear as the specified color, otherwise green + black
    symbol_size: int = 6
    line_width: int = 4
    tick_font: int = 14
    axis_title_font: int = 16
    plot_background_color: str = None  # Used to color the background of the 3D plot if desired
    camera_position: dict = None  # Determine if a specific camera value is used (See Plotly for Camera values)
    node_color: str = 'green'
    symbol: str = 'square'
    excluded_color: str = 'red'
    opaque_excluded_region: bool = False

    def __post_init__(self):
        self.fig = go.Figure()  # Initialize a figure to use in various functions
        if self.display_node_names:
            self.mode = "markers+text"
        else:
            self.mode = "markers"

        # Create the plotly layout figure:
        if self.display_background:
            background_color = self.plot_background_color
        else:
            background_color = 'rgba(0, 0, 0, 0)'

        if self.camera_position is None:
            camera_position = None
        else:
            camera_position = self.camera_position

        self.layout = go.Layout(
            scene=dict(
                xaxis_visible=self.display_axes,
                yaxis_visible=self.display_axes,
                zaxis_visible=self.display_axes,
                xaxis=dict(
                    title="X Axis (nm)",
                    tickfont=dict(size=self.tick_font, family='Helvetica', color='black'),
                    titlefont=dict(size=self.axis_title_font, family='Helvetica', color='black'),
                    showbackground=self.display_background,
                    showgrid=False,
                    backgroundcolor=background_color,
                    zerolinewidth=0.
                ),
                yaxis=dict(
                    title="Y Axis (nm)",
                    tickfont=dict(size=self.tick_font, family='Helvetica', color='black'),
                    titlefont=dict(size=self.axis_title_font, family='Helvetica', color='black'),
                    showbackground=self.display_background,
                    showgrid=False,
                    backgroundcolor=background_color,
                    zerolinewidth=0.
                ),
                zaxis=dict(
                    title="Z Axis (nm)",
                    tickfont=dict(size=self.tick_font, family='Helvetica', color='black'),
                    titlefont=dict(size=self.axis_title_font, family='Helvetica', color='black'),
                    showbackground=self.display_background,
                    showgrid=False,
                    backgroundcolor=background_color,
                    zerolinewidth=0.
                ),
                camera=camera_position
            ),
            plot_bgcolor="rgba(0, 0, 0, 0)",  # Set plot background color to white
            paper_bgcolor="rgba(0, 0, 0, 0)",  # Set paper background color to white
            margin=dict(l=0, r=0, b=0, t=0)
        )

        self.create_plot()  # Call the plot creation so the user can simply return or show the plot


    def plot_bounding_box(self) -> list:
        """ This method plots all bounding box edges """
        bounding_box_traces = []
        # Plotting Nodes:
        for node in self.bounding_box.bounding_graph.nodes():
            n = self.bounding_box.bounding_graph.nodes[node]
            x, y, z = n['x'], n['y'], n['z']
            node_trace = go.Scatter3d(x=[x], y=[y], z=[z], mode=self.mode, text=f'{node}',
                                      marker=dict(size=self.symbol_size, color='black', symbol='square'),
                                      textposition='top center', name='Bounding Region',
                                      legendgroup='Bounding Region', textfont=dict(size=12), showlegend=False)
            bounding_box_traces.append(node_trace)

        # Plotting Edges:
        for edge in self.bounding_box.bounding_graph.edges():
            n1, n2 = self.bounding_box.bounding_graph.nodes[edge[0]], self.bounding_box.bounding_graph.nodes[edge[1]]
            x1, y1, z1 = n1['x'], n1['y'], n1['z']
            x2, y2, z2 = n2['x'], n2['y'], n2['z']
            node_trace = go.Scatter3d(x=[x1, x2], y=[y1, y2], z=[z1, z2], mode='lines',
                                      line=dict(color='black', dash='dash', width=self.line_width),
                                      name='Bounding Region', legendgroup='Bounding Region', showlegend=False)
            bounding_box_traces.append(node_trace)
        return bounding_box_traces


    def plot_mesh_faces(self) -> list:
        """ This method creates all plotly traces for the mesh """
        all_face_traces = []
        for face in self.design.all_faces.values():
            if self.face_color is not None:
                trace = face.get_face_plot(color_type='manual', override_color=self.face_color)
            else:
                trace = face.get_face_plot(color_type='other')
            trace.legendgroup = 'dna'
            all_face_traces.append(trace)
        return all_face_traces


    def plot_edges_and_nodes(self) -> list:
        """ This method creates all traces containing the edges and nodes of the design """
        surface_mesh_traces = []
        if self.show_mesh_nodes:
            for node in self.design.design_graph.nodes():
                n = self.design.design_graph.nodes[node]
                x1, y1, z1 = n['x'], n['y'], n['z']
                if 'preserved_region' in n:
                    node_trace = go.Scatter3d(x=[x1], y=[y1], z=[z1], mode=self.mode, text=f'{node}',
                                              marker=dict(color=self.node_color, size=self.symbol_size, symbol=self.symbol),
                                              textposition='top center', name='Mesh Node',
                                              legendgroup='Graph Node', showlegend=False)
                else:
                    node_trace = go.Scatter3d(x=[x1], y=[y1], z=[z1], mode=self.mode, text=f'{node}',
                                              marker=dict(color='black', size=self.symbol_size, symbol=self.symbol),
                                              textposition='top center', name='Mesh Node',
                                              legendgroup='Graph Node', showlegend=False)
                surface_mesh_traces.append(node_trace)

        # Plot edges:
        for edge in self.design.design_graph.edges():
            n1, n2 = self.design.design_graph.nodes[edge[0]], self.design.design_graph.nodes[edge[1]]
            x1, y1, z1 = n1['x'], n1['y'], n1['z']
            x2, y2, z2 = n2['x'], n2['y'], n2['z']
            if self.constant_edge_color is not None:
                edge_trace = go.Scatter3d(x=[x1, x2], y=[y1, y2], z=[z1, z2], mode='lines',
                                          line=dict(color=self.constant_edge_color, width=self.line_width),
                                          name='Mesh Edge', legendgroup='Graph Edge', showlegend=False)

            # If both nodes are preserved_regions then we make the line green by default:
            elif edge in self.design.terminal_edges or (edge[1], edge[0]) in self.design.terminal_edges:
                edge_trace = go.Scatter3d(x=[x1, x2], y=[y1, y2], z=[z1, z2], mode='lines',
                                          line=dict(color=self.node_color, width=self.line_width), showlegend=False,
                                          name='Mesh Edge', legendgroup='Graph Edge')

            # Otherwise we use black for the line color:
            else:
                edge_trace = go.Scatter3d(x=[x1, x2], y=[y1, y2], z=[z1, z2], mode='lines',
                                          line=dict(color='black', width=self.line_width), showlegend=False,
                                          name='Mesh Edge', legendgroup='Graph Edge')

            surface_mesh_traces.append(edge_trace)


        return surface_mesh_traces


    def plot_excluded_regions(self) -> list:
        """ This method plots the excluded regions """
        excluded_face_traces = []
        for excluded_region in self.design.excluded:
            for face in excluded_region.faces.values():
                new_face_trace = face.get_face_plot(color_type='manual', override_opacity=self.opaque_excluded_region, override_color=self.excluded_color)
                new_face_trace.legendgroup = 'excluded'
                excluded_face_traces.append(new_face_trace)
        return excluded_face_traces


    def create_plot(self):
        """ This method is initialized by post_init and calls the above methods. Creates the plot. """
        all_traces = []
        if self.show_bounding_box:
            bounding_box_traces = self.plot_bounding_box()
            all_traces.extend(bounding_box_traces)

        if self.show_excluded_regions:
            excluded_traces = self.plot_excluded_regions()
            all_traces.extend(excluded_traces)

        # Grab the other traces:
        face_traces = self.plot_mesh_faces()
        all_traces.extend(face_traces)
        edge_n_node_traces = self.plot_edges_and_nodes()
        all_traces.extend(edge_n_node_traces)

        # Now we add to the figure:
        self.fig.add_traces(data=all_traces)
        self.fig.layout = self.layout  # Set the layout to the figure layout


    def return_figure(self):
        """ Method allowing user to use / manipulate the plotly figure of a design """
        return self.fig

    def show_figure(self):
        """ Function that shows the figure in the web browser """
        self.fig.show()


@dataclass
class CylindricalRepresentation(object):
    """
    Other than a mesh representation, this cylindrical representation is a bit more meaningful for DNA origami
    representations, since realistically the edge in the mesh dictates where material is at in the generated design.

    Parameters
    ------------
    design : Design we are visualizing of the PolyhedralSpace data object (see design_spaces.polyhedral_design_space).
    bounding_box : The bounding box defined, can be accessed by the design.bounding_box

    BOOLEAN options that turn on / off various visualizations:
    show_bounding_box : Determines if bounding box is shown, used for debugging (default: TRUE)
    show_excluded_regions : Determines if the excluded regions are shown or not (default: TRUE)
    opaque_excluded_region : Determines if excluded regions are opaque or not (default: FALSE / not opaque)
    display_background : Determines if a background color is to be used (default: False / blank background)
    display_axes : Determines if the axes and ticker labels are shown (default: True)
    show_box_bounds : Determines if custom annotations are displayed showing bounding box dimensions (default: FALSE)

    DISPLAY options that control colors / font sizes:
    cylinder_color : String color of the cylinders mapped by either plotly color or 'rgb(R, G, B)' (default: darkgray)
    cylinder_diameter : Floating point value of the cylinder diameter to use (default : 2.25) (BDNA diameter)
    tick_font : Integer point-size of ticker font (default: 14)
    axis_title_font : Integer point-size of axis titles
    camera_position : Camera dictionary (see plotly) for holding the camera angle constant
    excluded_color : String color of excluded regions by either plotly color or 'rgb(R, G, B)' (default: red)
    plot_background_color : String color of background defined by either plotly color or 'rgb(R, G, B)' (default: None)
    """
    design: PolyhedralSpace
    bounding_box: Union[TriclinicBox, MonoclinicBox, OrthorhombicBox, TetragonalBox,
                        RhombohedralBox, HexagonalBox, CubicBox]
    show_bounding_box: bool = True
    show_box_bounds: bool = False
    show_excluded_regions: bool = True
    display_background: bool = False
    display_axes: bool = True  # Determine if axis titles and ticker labels are shown
    cylinder_color: str = 'darkgray'  # Override potential face color of mesh if desired
    excluded_color: str = 'red'
    tick_font: int = 18
    axis_title_font: int = 22
    plot_background_color: str = None  # Used to color the background of the 3D plot if desired
    camera_position: dict = None  # Determine if a specific camera value is used (See Plotly for Camera values)
    cylinder_diameter: float = 2.25  # Diameter of DNA to use to represent as cylinder
    opaque_excluded_region: bool = False

    def __post_init__(self):
        self.fig = go.Figure()  # Initialize a figure to use in various functions

        # Create the plotly layout figure:
        if self.display_background:
            background_color = self.plot_background_color
        else:
            background_color = 'rgba(0, 0, 0, 0)'

        if self.camera_position is None:
            camera_position = None
        else:
            camera_position = self.camera_position

        self.layout = go.Layout(
            scene=dict(
                xaxis_visible=self.display_axes,
                yaxis_visible=self.display_axes,
                zaxis_visible=self.display_axes,
                xaxis=dict(
                    title="X Axis (nm)",
                    tickfont=dict(size=self.tick_font, family='Helvetica', color='black'),
                    titlefont=dict(size=self.axis_title_font, family='Helvetica', color='black'),
                    showbackground=self.display_background,
                    showgrid=False,
                    backgroundcolor=background_color,
                    zerolinewidth=0.
                ),
                yaxis=dict(
                    title="Y Axis (nm)",
                    tickfont=dict(size=self.tick_font, family='Helvetica', color='black'),
                    titlefont=dict(size=self.axis_title_font, family='Helvetica', color='black'),
                    showbackground=self.display_background,
                    showgrid=False,
                    backgroundcolor=background_color,
                    zerolinewidth=0.
                ),
                zaxis=dict(
                    title="Z Axis (nm)",
                    tickfont=dict(size=self.tick_font, family='Helvetica', color='black'),
                    titlefont=dict(size=self.axis_title_font, family='Helvetica', color='black'),
                    showbackground=self.display_background,
                    showgrid=False,
                    backgroundcolor=background_color,
                    zerolinewidth=0.
                ),
                camera=camera_position
            ),
            plot_bgcolor="rgba(0, 0, 0, 0)",  # Set plot background color to white
            paper_bgcolor="rgba(0, 0, 0, 0)",  # Set paper background color to white
            margin=dict(l=0, r=0, b=0, t=0)
        )

        self.create_plot()  # Call the plot creation so the user can simply return or show the plot


    def plot_bounding_box(self) -> list:
        """ This method plots all bounding box edges """
        bounding_box_traces = []

        # Plotting Edges:
        for edge in self.bounding_box.bounding_graph.edges():
            n1, n2 = self.bounding_box.bounding_graph.nodes[edge[0]], self.bounding_box.bounding_graph.nodes[edge[1]]
            x1, y1, z1 = n1['x'], n1['y'], n1['z']
            x2, y2, z2 = n2['x'], n2['y'], n2['z']
            node_trace = go.Scatter3d(x=[x1, x2], y=[y1, y2], z=[z1, z2], mode='lines',
                                      line=dict(color='black', dash='dash', width=6),
                                      name='Bounding Region', legendgroup='Bounding Region', showlegend=False)
            bounding_box_traces.append(node_trace)
        return bounding_box_traces


    def plot_edges_as_surfaces(self) -> list:
        """ This method plots all edges using a cylindrical approximation """
        surface_mesh_traces = []
        CYLINDER_RESOLUTION = 25  # Testing with the value and 25 seems to be fine computationally speaking
        r = self.cylinder_diameter / 2

        for edge in self.design.design_graph.edges():
            length, direction = length_and_direction_between_nodes(graph=self.design.design_graph,
                                                                   node1=edge[0], node2=edge[1])

            rotation_matrix = rotation_matrix_from_vectors(np.array([0, 0, 1]), direction)
            theta = np.linspace(0, 2 * np.pi, CYLINDER_RESOLUTION)
            z = np.linspace(0, length, CYLINDER_RESOLUTION)
            theta_grid, z_grid = np.meshgrid(theta, z)
            cylinder_surface_points = np.array([r * np.cos(theta_grid), r * np.sin(theta_grid), z_grid])
            if np.allclose(rotation_matrix, np.eye(3)):
                # Handle the case when vectors are almost parallel
                P1 = find_smaller_z_value(graph=self.design.design_graph, node1=edge[0], node2=edge[1])
            else:
                # Apply the rotation matrix to cylinder surface points
                cylinder_surface_points = np.array([r * np.cos(theta_grid), r * np.sin(theta_grid), z_grid])
                for i in range(CYLINDER_RESOLUTION):
                    for j in range(CYLINDER_RESOLUTION):
                        cylinder_surface_points[:, i, j] = np.dot(rotation_matrix, cylinder_surface_points[:, i, j])
                # Translate the points to the desired position in 3D space
                n1 = self.design.design_graph.nodes[edge[0]]
                P1 = np.array([n1['x'], n1['y'], n1['z']])

            cylinder_surface_points = cylinder_surface_points + P1.reshape(-1, 1, 1)
            x_grid, y_grid, z_grid = cylinder_surface_points[0], cylinder_surface_points[1], cylinder_surface_points[2]
            x_points = x_grid.reshape(-1)
            y_points = y_grid.reshape(-1)
            z_points = z_grid.reshape(-1)

            # Create vertices for the cylinder's surface
            faces, vertices = [], []
            for i in range(CYLINDER_RESOLUTION):
                for j in range(CYLINDER_RESOLUTION):
                    x = x_points[i * CYLINDER_RESOLUTION + j]
                    y = y_points[i * CYLINDER_RESOLUTION + j]
                    z = z_points[i * CYLINDER_RESOLUTION + j]
                    vertices.append([x, y, z])

            # Create faces by defining the vertex indices
            for i in range(CYLINDER_RESOLUTION - 1):
                for j in range(CYLINDER_RESOLUTION - 1):
                    v1 = i * CYLINDER_RESOLUTION + j
                    v2 = v1 + 1
                    v3 = (i + 1) * CYLINDER_RESOLUTION + j
                    v4 = v3 + 1
                    faces.extend([(v1, v2, v3), (v2, v4, v3)])

            new_trace = go.Mesh3d(x=x_points, y=y_points, z=z_points, i=[face[0] for face in faces],
                                  j=[face[1] for face in faces], k=[face[2] for face in faces],
                                  color=self.cylinder_color, showlegend=False, legendgroup='dna')
            surface_mesh_traces.append(new_trace)

        return surface_mesh_traces


    def plot_excluded_regions(self) -> list:
        """ This method plots the excluded regions """
        excluded_face_traces = []
        for excluded_region in self.design.excluded:
            for face in excluded_region.faces.values():
                new_face_trace = face.get_face_plot(color_type='manual', override_opacity=self.opaque_excluded_region,
                                                    override_color=self.excluded_color)
                new_face_trace.legendgroup = 'excluded'
                excluded_face_traces.append(new_face_trace)
        return excluded_face_traces


    def create_plot(self):
        """ This method is initialized by post_init and calls the above methods. Creates the plot. """
        all_traces = []
        if self.show_bounding_box:
            bounding_box_traces = self.plot_bounding_box()
            all_traces.extend(bounding_box_traces)

        if self.show_box_bounds:
            # This add the "Unit Cell" Parameters to the plot
            cone_traces, trace_annotations, text_annotation = self.bounding_box.return_box_traces(text_size=self.axis_title_font+4)
            # Add text annotations
            all_traces.extend(cone_traces)
            all_traces.extend(trace_annotations)
            if text_annotation is not None:
                self.layout.scene.annotations = [text_annotation]


        if self.show_excluded_regions:
            excluded_traces = self.plot_excluded_regions()
            all_traces.extend(excluded_traces)

        # Grab the other traces:
        edge_n_node_traces = self.plot_edges_as_surfaces()
        all_traces.extend(edge_n_node_traces)

        # Now we add to the figure:
        self.fig.add_traces(data=all_traces)
        self.fig.layout = self.layout  # Set the layout to the figure layout


    def return_figure(self):
        """ Method allowing user to use / manipulate the plotly figure of a design """
        return self.fig

    def show_figure(self):
        """ Function that shows the figure in the web browser """
        self.fig.show()
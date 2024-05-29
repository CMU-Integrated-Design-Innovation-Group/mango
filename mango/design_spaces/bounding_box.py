"""
A.J. Vetturini
IDIG and MMBL
Carnegie Mellon University

This module defines the bounding box surrounding the polyhedral design space, it has many classes for the variety of
box-types that the bounding_box can presume.
"""
from dataclasses import dataclass
from mango.utils.mango_math import *
import plotly.graph_objs as go

# Function call use in boxes:
def create_parallelepiped(cell_type: str, a: float, b: float, c: float, alpha: float, beta: float, gamma: float):
    return Parallelepiped(cell_type=cell_type, a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma)


@dataclass
class TriclinicBox(object):
    """
    This is a bounding box exhibiting Triclinic Conditions by default. It is the most "flexible" space as every
    parameter of the parallelepiped can be constrained.

    Parameters
    ------------
    a : Magnitude of parallelepiped aligned in x
    b : Magnitude of parallelepiped aligned in y
    c: Magnitude of parallelepiped aligned in z
    alpha : Angle of parallelepiped between b and c
    beta : Angle of parallelepiped between c and a
    gamma: Angle of parallelepiped between b and a
    """
    a: float  # Aligns in X
    b: float   # Aligns in Y
    c: float  # Aligns in Z

    # Given in units of degrees:
    alpha: float  # Angle between b and c
    beta: float  # Angle between c and a
    gamma: float  # Angle between b and a

    def __post_init__(self):
        if not self.edges_unique():
            raise Exception('Invalid edge values for a triclinic box, a b and c must be unique.')
        if not self.angles_unique():
            raise Exception('Invalid edge values for a triclinic box, alpha beta and gamma must be unique.')

        for i in [self.a, self.b, self.c]:
            if not positive_value(i):
                raise Exception('A specified value for a, b, or c is less than zero and therefore invalid. Please'
                                ' check your entered values.')

        # Create box as parallelepiped
        self.shape = create_parallelepiped(cell_type='triclinic', a=self.a, b=self.b, c=self.c, alpha=self.alpha,
                                           beta=self.beta, gamma=self.gamma)
        self.bounding_graph = self.shape.bounding_graph

        # Map of values:
        self.param_map = {
            'a': ['a'],
            'b': ['b'],
            'c': ['c'],
            'alpha': ['alpha'],
            'beta': ['beta'],
            'gamma': ['gamma']
        }

    def edges_unique(self) -> bool:
        """
        Validates if any one edge is the same value
        """
        if self.a == self.b or self.a == self.c or self.b == self.c:
            return False
        return True

    def angles_unique(self) -> bool:
        """
        Validates if any one angle is the same value
        """
        if self.alpha == self.beta or self.alpha == self.gamma or self.beta == self.gamma:
            return False
        return True

    def update_shape(self) -> None:
        """
        Function that simple updated the shape and bounding_graphs based on the updated values of a, b, c
        """
        self.shape = create_parallelepiped(cell_type='triclinic', a=self.a, b=self.b, c=self.c, alpha=self.alpha,
                                           beta=self.beta, gamma=self.gamma)
        self.bounding_graph = self.shape.bounding_graph  # We also pass back in / updated the bounding graph

    def return_box_traces(self, text_size: float) -> tuple:
        """
        This function returns traces displaying the a b c and alpha beta gamma values via plotly

        :param text_size: Desired font size in units of pt
        :return: Various traces that are assigned to a plotly figure
        """
        # Next create the traces for the box given Triclinic box:
        cone_traces = [
            go.Cone(x=[10], y=[1], z=[1], u=[1], v=[0], w=[0], sizeref=12,
                    colorscale=[[0, 'rgb(187, 85, 102)'], [1, 'rgb(187, 85, 102)']], showscale=False, showlegend=False),
            go.Cone(x=[1], y=[10], z=[1], u=[0], v=[1], w=[0], sizeref=12,
                    colorscale=[[0, 'rgb(221, 170, 51)'], [1, 'rgb(221, 170, 51)']], showscale=False, showlegend=False),
            go.Cone(x=[1], y=[1], z=[10], u=[0], v=[0], w=[1], sizeref=12,
                    colorscale=[[0, 'rgb(0, 68, 136)'], [1, 'rgb(0, 68, 136)']], showscale=False, showlegend=False),
        ]

        annotation_text = f'\u03B1<sub style="color:rgb(221, 170, 51);">b</sub><sub style="color:rgb(0, 68, 136);">c</sub> : {self.alpha:.2f}\u00b0<br>' \
                          f'\u03B2<sub style="color:rgb(187, 85, 102);">a</sub><sub style="color:rgb(0, 68, 136);">c</sub>: {self.beta:.2f}\u00b0<br>' \
                          f'\u03B3<sub style="color:rgb(187, 85, 102);">a</sub><sub style="color:rgb(221, 170, 51);">b</sub> : {self.gamma:.2f}\u00b0'

        text_annotation = dict(
            x=0.1,  # Adjust x position as needed
            y=0.9,  # Adjust y position as needed
            z=self.c,
            text=annotation_text,
            showarrow=False,
            font=dict(
                family="Helvetica",
                size=text_size,
                color="black"
            )
        )

        trace_annotations = [
            # Add values of a, b, c
            go.Scatter3d(x=[20], y=[3], z=[3], mode='text', text=[f'a = {round(self.a, 2)}'],
                         textfont=dict(size=text_size, color='rgb(187, 85, 102)'),
                         showlegend=False, textposition="top center", hoverinfo='none'),
            go.Scatter3d(x=[3], y=[20], z=[3], mode='text', text=[f'b = {round(self.b, 2)}'],
                         textfont=dict(size=text_size, color='rgb(221, 170, 51)'),
                         showlegend=False, textposition="top center", hoverinfo='none'),
            go.Scatter3d(x=[3], y=[3], z=[20], mode='text', text=[f'c = {round(self.c, 2)}'],
                         textfont=dict(size=text_size, color='rgb(0, 68, 136)'),
                         showlegend=False, textposition="top center", hoverinfo='none')
        ]
        return cone_traces, trace_annotations, text_annotation



@dataclass
class MonoclinicBox(object):
    """
    This is a bounding box exhibiting Monoclinic Conditions by default. It's edge lengths can be different along with
    an angle specified between c and a, but alpha and gamma are held to 90 degrees by default

    Parameters
    ------------
    a : Magnitude of parallelepiped aligned in x
    b : Magnitude of parallelepiped aligned in y
    c: Magnitude of parallelepiped aligned in z
    beta : Angle of parallelepiped between c and a
    """
    a: float
    b: float
    c: float

    # Given in units of degrees:
    beta: float

    def __post_init__(self):
        self.alpha = 90
        self.gamma = 90
        # Map of values:
        self.param_map = {
            'a': ['a'],
            'b': ['b'],
            'c': ['c'],
            'beta': ['beta'],
        }

        if not self.edges_unique():
            raise Exception('Invalid edge values for a monoclinic box, a b and c must be unique.')

        for i in [self.a, self.b, self.c]:
            if not positive_value(i):
                raise Exception('A specified value for a, b, or c is less than zero and therefore invalid. Please'
                                ' check your entered values.')

        # Create box as parallelepiped
        self.shape = create_parallelepiped(cell_type='monoclinic', a=self.a, b=self.b, c=self.c, alpha=self.alpha,
                                           beta=self.beta, gamma=self.gamma)
        self.bounding_graph = self.shape.bounding_graph

    def edges_unique(self) -> bool:
        """
        Validates if any one edge is the same value
        """
        if self.a == self.b or self.a == self.c or self.b == self.c:
            return False
        return True

    def update_shape(self):
        """
        Function that simple updated the shape and bounding_graphs based on the updated values
        """
        self.shape = create_parallelepiped(cell_type='monoclinic', a=self.a, b=self.b, c=self.c, alpha=self.alpha,
                                           beta=self.beta, gamma=self.gamma)
        self.bounding_graph = self.shape.bounding_graph  # We also pass back in / updated the bounding graph

    def return_box_traces(self, text_size: float) -> tuple:
        """
        This function returns traces displaying the a b c and alpha beta gamma values via plotly

        :param text_size: Desired font size in units of pt
        :return: Various traces that are assigned to a plotly figure
        """
        # Next create the traces for the box given Triclinic box:
        cone_traces = [
            go.Cone(x=[10], y=[1], z=[1], u=[1], v=[0], w=[0], sizeref=12,
                    colorscale=[[0, 'rgb(187, 85, 102)'], [1, 'rgb(187, 85, 102)']], showscale=False, showlegend=False),
            go.Cone(x=[1], y=[10], z=[1], u=[0], v=[1], w=[0], sizeref=12,
                    colorscale=[[0, 'rgb(221, 170, 51)'], [1, 'rgb(221, 170, 51)']], showscale=False, showlegend=False),
            go.Cone(x=[1], y=[1], z=[10], u=[0], v=[0], w=[1], sizeref=12,
                    colorscale=[[0, 'rgb(0, 68, 136)'], [1, 'rgb(0, 68, 136)']], showscale=False, showlegend=False),
        ]

        annotation_text = f'\u03B2<sub style="color:rgb(187, 85, 102);">a</sub><sub style="color:rgb(0, 68, 136);">c</sub>: {self.beta:.2f}\u00b0<br>' \

        text_annotation = dict(
            x=0.1,  # Adjust x position as needed
            y=0.9,  # Adjust y position as needed
            z=self.c,
            text=annotation_text,
            showarrow=False,
            font=dict(
                family="Helvetica",
                size=text_size,
                color="black"
            )
        )

        trace_annotations = [
            # Add values of a, b, c
            go.Scatter3d(x=[20], y=[3], z=[3], mode='text', text=[f'a = {round(self.a, 2)}'],
                         textfont=dict(size=text_size, color='rgb(187, 85, 102)'),
                         showlegend=False, textposition="top center", hoverinfo='none'),
            go.Scatter3d(x=[3], y=[20], z=[3], mode='text', text=[f'b = {round(self.b, 2)}'],
                         textfont=dict(size=text_size, color='rgb(221, 170, 51)'),
                         showlegend=False, textposition="top center", hoverinfo='none'),
            go.Scatter3d(x=[3], y=[3], z=[20], mode='text', text=[f'c = {round(self.c, 2)}'],
                         textfont=dict(size=text_size, color='rgb(0, 68, 136)'),
                         showlegend=False, textposition="top center", hoverinfo='none')
        ]
        return cone_traces, trace_annotations, text_annotation


@dataclass
class OrthorhombicBox(object):
    """
    This is a bounding box exhibiting Orthorhombic Conditions by default. It's edge lengths can be different but the
    defining angles of the box are all held perpendicular (90 degrees).

    Parameters
    ------------
    a : Magnitude of parallelepiped aligned in x
    b : Magnitude of parallelepiped aligned in y
    c: Magnitude of parallelepiped aligned in z
    """
    a: float
    b: float
    c: float

    def __post_init__(self):
        self.alpha = 90
        self.beta = 90
        self.gamma = 90
        self.param_map = {
            'a': ['a'],
            'b': ['b'],
            'c': ['c'],
        }

        if not self.edges_unique():
            raise Exception('Invalid edge values for a orthorhombic box, a b and c must be unique.')

        for i in [self.a, self.b, self.c]:
            if not positive_value(i):
                raise Exception('A specified value for a, b, or c is less than zero and therefore invalid. Please'
                                ' check your entered values.')

        # Create box as parallelepiped
        self.shape = create_parallelepiped(cell_type='orthorhombic', a=self.a, b=self.b, c=self.c, alpha=self.alpha,
                                           beta=self.beta, gamma=self.gamma)
        self.bounding_graph = self.shape.bounding_graph

    def edges_unique(self) -> bool:
        if self.a == self.b or self.a == self.c or self.b == self.c:
            return False
        return True

    def update_shape(self):
        """
        Function that simple updated the shape and bounding_graphs based on the updated values
        """
        self.shape = create_parallelepiped(cell_type='orthorhombic', a=self.a, b=self.b, c=self.c, alpha=self.alpha,
                                           beta=self.beta, gamma=self.gamma)
        self.bounding_graph = self.shape.bounding_graph  # We also pass back in / updated the bounding graph

    def return_box_traces(self, text_size: float) -> tuple:
        """
        This function returns traces displaying the a b c and alpha beta gamma values via plotly

        :param text_size: Desired font size in units of pt
        :return: Various traces that are assigned to a plotly figure
        """
        # Next create the traces for the box given Triclinic box:
        cone_traces = [
            go.Cone(x=[10], y=[1], z=[1], u=[1], v=[0], w=[0], sizeref=12,
                    colorscale=[[0, 'rgb(187, 85, 102)'], [1, 'rgb(187, 85, 102)']], showscale=False, showlegend=False),
            go.Cone(x=[1], y=[10], z=[1], u=[0], v=[1], w=[0], sizeref=12,
                    colorscale=[[0, 'rgb(221, 170, 51)'], [1, 'rgb(221, 170, 51)']], showscale=False, showlegend=False),
            go.Cone(x=[1], y=[1], z=[10], u=[0], v=[0], w=[1], sizeref=12,
                    colorscale=[[0, 'rgb(0, 68, 136)'], [1, 'rgb(0, 68, 136)']], showscale=False, showlegend=False),
        ]

        trace_annotations = [
            go.Scatter3d(x=[20], y=[3], z=[3], mode='text', text=[f'a = {round(self.a, 2)}'],
                         textfont=dict(size=text_size, color='rgb(187, 85, 102)'),
                         showlegend=False, textposition="top center", hoverinfo='none'),
            go.Scatter3d(x=[3], y=[20], z=[3], mode='text', text=[f'b = {round(self.b, 2)}'],
                         textfont=dict(size=text_size, color='rgb(221, 170, 51)'),
                         showlegend=False, textposition="top center", hoverinfo='none'),
            go.Scatter3d(x=[3], y=[3], z=[20], mode='text', text=[f'c = {round(self.c, 2)}'],
                         textfont=dict(size=text_size, color='rgb(0, 68, 136)'),
                         showlegend=False, textposition="top center", hoverinfo='none')
        ]
        return cone_traces, trace_annotations, None

@dataclass
class TetragonalBox(object):
    """
    This is a bounding box exhibiting Tetragonal Conditions by default. Here, a and b are held to the same value and c
    defines a height. All angles are held perpendicular

    Parameters
    ------------
    a : Magnitude of parallelepiped aligned in x, y
    c: Magnitude of parallelepiped aligned in z
    """
    a: float
    c: float

    def __post_init__(self):
        self.alpha = 90
        self.beta = 90
        self.gamma = 90
        self.b = self.a

        self.param_map = {
            'a': ['a', 'b'],
            'c': ['c'],
        }

        if self.a == self.c:
            raise Exception('Invalid edge values for a tetragonal box, a and c must be unique.')

        for i in [self.a, self.b, self.c]:
            if not positive_value(i):
                raise Exception('A specified value for a, b, or c is less than zero and therefore invalid. Please'
                                ' check your entered values.')

        # Create box as parallelepiped
        self.shape = create_parallelepiped(cell_type='tetragonal', a=self.a, b=self.b, c=self.c, alpha=self.alpha,
                                           beta=self.beta, gamma=self.gamma)
        self.bounding_graph = self.shape.bounding_graph

    def update_shape(self):
        """
        Function that simple updated the shape and bounding_graphs based on the updated values
        """
        self.shape = create_parallelepiped(cell_type='tetragonal', a=self.a, b=self.b, c=self.c, alpha=self.alpha,
                                           beta=self.beta, gamma=self.gamma)
        self.bounding_graph = self.shape.bounding_graph  # We also pass back in / updated the bounding graph

    def return_box_traces(self, text_size: float) -> tuple:
        """
        This function returns traces displaying the a b c and alpha beta gamma values via plotly

        :param text_size: Desired font size in units of pt
        :return: Various traces that are assigned to a plotly figure
        """
        # Next create the traces for the box given Triclinic box:
        cone_traces = [
            go.Cone(x=[10], y=[1], z=[1], u=[1], v=[0], w=[0], sizeref=12,
                    colorscale=[[0, 'rgb(187, 85, 102)'], [1, 'rgb(187, 85, 102)']], showscale=False, showlegend=False),
            go.Cone(x=[1], y=[1], z=[10], u=[0], v=[0], w=[1], sizeref=12,
                    colorscale=[[0, 'rgb(0, 68, 136)'], [1, 'rgb(0, 68, 136)']], showscale=False, showlegend=False),
        ]

        trace_annotations = [
            # Add values of a, c
            go.Scatter3d(x=[20], y=[3], z=[3], mode='text', text=[f'a = {round(self.a, 2)}'],
                         textfont=dict(size=text_size, color='rgb(187, 85, 102)'),
                         hoverlabel=dict(bgcolor='rgb(128, 128, 128)'),
                         showlegend=False, textposition="top center", hoverinfo='none'),
            go.Scatter3d(x=[3], y=[3], z=[20], mode='text', text=[f'c = {round(self.c, 2)}'],
                         textfont=dict(size=text_size, color='rgb(0, 68, 136)'),
                         showlegend=False, textposition="top center", hoverinfo='none')
        ]
        return cone_traces, trace_annotations, None


@dataclass
class RhombohedralBox(object):
    """
    This is a bounding box exhibiting Rhombohedral Conditions by default. Here, all edges and angles are held to the
    same values, respectively.

    Parameters
    ------------
    a : Magnitude of all parallelepiped edges
    alpha : Angle of parallelepiped between all edges
    """
    a: float
    alpha: float

    def __post_init__(self):
        self.beta = self.alpha
        self.gamma = self.alpha
        self.b = self.a
        self.c = self.a

        self.param_map = {
            'a': ['a', 'b', 'c'],
            'alpha': ['alpha', 'beta', 'gamma'],
        }

        for i in [self.a, self.b, self.c]:
            if not positive_value(i):
                raise Exception('A specified value for a, b, or c is less than zero and therefore invalid. Please'
                                ' check your entered values.')

        # Create box as parallelepiped
        self.shape = create_parallelepiped(cell_type='rhombohedral', a=self.a, b=self.b, c=self.c, alpha=self.alpha,
                                           beta=self.beta, gamma=self.gamma)
        self.bounding_graph = self.shape.bounding_graph

    def update_shape(self):
        """
        Function that simple updated the shape and bounding_graphs based on the updated values
        """
        self.shape = create_parallelepiped(cell_type='rhombohedral', a=self.a, b=self.b, c=self.c, alpha=self.alpha,
                                           beta=self.beta, gamma=self.gamma)
        self.bounding_graph = self.shape.bounding_graph  # We also pass back in / updated the bounding graph

    def return_box_traces(self, text_size: float) -> tuple:
        """
        This function returns traces displaying the a b c and alpha beta gamma values via plotly

        :param text_size: Desired font size in units of pt
        :return: Various traces that are assigned to a plotly figure
        """
        # Next create the traces for the box given Triclinic box:
        cone_traces = [
            go.Cone(x=[10], y=[1], z=[1], u=[1], v=[0], w=[0], sizeref=12,
                    colorscale=[[0, 'rgb(187, 85, 102)'], [1, 'rgb(187, 85, 102)']], showscale=False, showlegend=False),
        ]

        annotation_text = f'\u03B1<sub style="color:rgb(221, 170, 51);">b</sub><sub style="color:rgb(0, 68, 136);">c</sub> : {self.alpha:.2f}\u00b0<br>'

        text_annotation = dict(
            x=0.1,  # Adjust x position as needed
            y=0.9,  # Adjust y position as needed
            z=self.c,
            text=annotation_text,
            showarrow=False,
            font=dict(
                family="Helvetica",
                size=text_size,
                color="black"
            )
        )

        trace_annotations = [
            # Add values of a, b, c
            go.Scatter3d(x=[20], y=[3], z=[3], mode='text', text=[f'a = {round(self.a, 2)}'],
                         textfont=dict(size=text_size, color='rgb(187, 85, 102)'),
                         showlegend=False, textposition="top center", hoverinfo='none'),
        ]
        return cone_traces, trace_annotations, text_annotation

@dataclass
class HexagonalBox(object):
    """
    This is a bounding box exhibiting Hexagonal Conditions by default. Here the edges in X and Y are held the same
    and a value of the height is specified. Further, the angles alpha and beta are held to 60 whereas gamma is 120.

    Parameters
    ------------
    a : Magnitude of parallelepiped aligned in x, y
    c: Magnitude of parallelepiped aligned in z
    """
    a: float
    c: float

    def __post_init__(self):
        self.gamma = 120
        self.alpha = 60
        self.beta = 60
        self.b = self.a

        self.param_map = {
            'a': ['a', 'b'],
            'c': ['c'],
        }

        for i in [self.a, self.b, self.c]:
            if not positive_value(i):
                raise Exception('A specified value for a, b, or c is less than zero and therefore invalid. Please'
                                ' check your entered values.')

        # Create box as parallelepiped
        self.shape = create_parallelepiped(cell_type='hexagonal', a=self.a, b=self.b, c=self.c, alpha=self.alpha, beta=self.beta, gamma=self.gamma)
        self.bounding_graph = self.shape.bounding_graph

    def update_shape(self):
        """
        Function that simple updated the shape and bounding_graphs based on the updated values
        """
        self.shape = create_parallelepiped(cell_type='hexagonal', a=self.a, b=self.b, c=self.c, alpha=self.alpha,
                                           beta=self.beta, gamma=self.gamma)
        self.bounding_graph = self.shape.bounding_graph  # We also pass back in / updated the bounding graph

    def return_box_traces(self, text_size: float) -> tuple:
        """
        This function returns traces displaying the a b c and alpha beta gamma values via plotly

        :param text_size: Desired font size in units of pt
        :return: Various traces that are assigned to a plotly figure
        """
        # Next create the traces for the box given Triclinic box:
        cone_traces = [
            go.Cone(x=[10], y=[1], z=[1], u=[1], v=[0], w=[0], sizeref=12,
                    colorscale=[[0, 'rgb(187, 85, 102)'], [1, 'rgb(187, 85, 102)']], showscale=False, showlegend=False),
            go.Cone(x=[1], y=[1], z=[10], u=[0], v=[0], w=[1], sizeref=12,
                    colorscale=[[0, 'rgb(0, 68, 136)'], [1, 'rgb(0, 68, 136)']], showscale=False, showlegend=False),
        ]

        annotation_text = f'\u03B1<sub style="color:rgb(221, 170, 51);">b</sub><sub style="color:rgb(0, 68, 136);">c</sub> : {self.alpha:.2f}\u00b0<br>' \
                          f'\u03B2<sub style="color:rgb(187, 85, 102);">a</sub><sub style="color:rgb(0, 68, 136);">c</sub>: {self.beta:.2f}\u00b0<br>' \
                          f'\u03B3<sub style="color:rgb(187, 85, 102);">a</sub><sub style="color:rgb(221, 170, 51);">b</sub> : {self.gamma:.2f}\u00b0'

        text_annotation = dict(
            x=0.1,  # Adjust x position as needed
            y=0.9,  # Adjust y position as needed
            z=self.c,
            text=annotation_text,
            showarrow=False,
            font=dict(
                family="Helvetica",
                size=text_size,
                color="black"
            )
        )

        trace_annotations = [
            # Add values of a, b, c
            go.Scatter3d(x=[20], y=[3], z=[3], mode='text', text=[f'a = {round(self.a, 2)}'],
                         textfont=dict(size=text_size, color='rgb(187, 85, 102)'),
                         showlegend=False, textposition="top center", hoverinfo='none'),
            go.Scatter3d(x=[3], y=[3], z=[20], mode='text', text=[f'c = {round(self.c, 2)}'],
                         textfont=dict(size=text_size, color='rgb(0, 68, 136)'),
                         showlegend=False, textposition="top center", hoverinfo='none')

        ]
        return cone_traces, trace_annotations, text_annotation


@dataclass
class CubicBox(object):  # Ie isometric
    """
    This is a bounding box exhibiting Isometric (Cubic) Conditions by default. All edges are held to the same length and
    all angles are held to 90 degrees.

    Parameters
    ------------
    a : Magnitude of all parallelepiped edges
    """
    a: float

    def __post_init__(self):
        self.gamma = 90
        self.alpha = 90
        self.beta = 90
        self.b = self.a
        self.c = self.a
        self.param_map = {
            'a': ['a', 'b', 'c'],
        }

        for i in [self.a, self.b, self.c]:
            if not positive_value(i):
                raise Exception('A specified value for a, b, or c is less than zero and therefore invalid. Please'
                                ' check your entered values.')

        # Create box as parallelepiped
        self.shape = create_parallelepiped(cell_type='cubic', a=self.a, b=self.b, c=self.c, alpha=self.alpha,
                                           beta=self.beta, gamma=self.gamma)
        self.bounding_graph = self.shape.bounding_graph

    def update_shape(self):
        """
        Function that simple updated the shape and bounding_graphs based on the updated values
        """
        self.shape = create_parallelepiped(cell_type='cubic', a=self.a, b=self.b, c=self.c, alpha=self.alpha,
                                           beta=self.beta, gamma=self.gamma)
        self.bounding_graph = self.shape.bounding_graph  # We also pass back in / updated the bounding graph

    def return_box_traces(self, text_size: float) -> tuple:
        """
        This function returns traces displaying the a b c and alpha beta gamma values via plotly

        :param text_size: Desired font size in units of pt
        :return: Various traces that are assigned to a plotly figure
        """
        # Next create the traces for the box given Triclinic box:
        cone_traces = [
            go.Cone(x=[10], y=[1], z=[1], u=[1], v=[0], w=[0], sizeref=12,
                    colorscale=[[0, 'rgb(187, 85, 102)'], [1, 'rgb(187, 85, 102)']], showscale=False, showlegend=False),
        ]

        trace_annotations = [
            # Add values of a, b, c
            go.Scatter3d(x=[20], y=[3], z=[3], mode='text', text=[f'a = {round(self.a, 2)}'],
                         textfont=dict(size=text_size, color='rgb(187, 85, 102)'),
                         showlegend=False, textposition="top center", hoverinfo='none'),
        ]
        return cone_traces, trace_annotations, None



@dataclass
class Parallelepiped(object):
    """
    This is a high level class that is used to initialize the various default conditions.

    Parameters
    ------------
    a : Magnitude of parallelepiped aligned in x
    b : Magnitude of parallelepiped aligned in y
    c: Magnitude of parallelepiped aligned in z
    alpha : Angle of parallelepiped between b and c
    beta : Angle of parallelepiped between c and a
    gamma: Angle of parallelepiped between b and a
    cell_type : What system is being shown (used in calculating volume)
    angles_in_degrees : Default TRUE, If false: will presume alpha beta gamma are in radians
    """
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float
    cell_type: str
    angles_in_degrees: bool = True

    def __post_init__(self):
        # Due to computational geometry, we want to ensure the faces of the Parallelepiped are constant so we declare
        # a constant rotation here going CCW at each face:
        self.faces = [(0, 2, 6, 1), (3, 5, 7, 4), (0, 2, 5, 3), (1, 4, 7, 6), (2, 6, 7, 5), (0, 3, 4, 1)]

        # Update angles to be in rads:
        if self.angles_in_degrees:
            self.alpha = np.deg2rad(self.alpha)
            self.beta = np.deg2rad(self.beta)
            self.gamma = np.deg2rad(self.gamma)

        self.bounding_graph = nx.empty_graph()

        # Calculate vertices and edges:
        vertices = self.calculate_vertices_and_edges()
        nodes = [(f'b{i}', {'bounding': True, 'x': v[0], 'y': v[1], 'z': v[2]}) for i, v in enumerate(vertices, 1)]

        self.bounding_graph.add_nodes_from(nodes)
        edges = [('b1', 'b2'), ('b1', 'b3'), ('b1', 'b4'),
                 ('b2', 'b5'), ('b2', 'b7'), ('b3', 'b7'),
                 ('b3', 'b6'), ('b4', 'b6'), ('b4', 'b5'),
                 ('b5', 'b8'), ('b6', 'b8'), ('b7', 'b8')]
        self.bounding_graph.add_edges_from(edges)
        self.update_face_equations()
        self.volume = 0.
        self.update_volume()
        self.midpoints = self.get_midpoint_of_each_face(vertices=vertices)


    def get_midpoint_of_each_face(self, vertices: list) -> list:
        """
        This function will determine the mid-point of each face. We can then add grammars to "move" these face-based
        points if desired, for now we will hold this constant.

        :param vertices: List of vertex node numbers we can access to pull correct values from graph
        :return: List of midpoints in any face
        """

        def midpoint_of_face(verts):
            # Find the midpoint of a face given its vertices
            x_sum = sum(vertex[0] for vertex in verts)
            y_sum = sum(vertex[1] for vertex in verts)
            z_sum = sum(vertex[2] for vertex in verts)

            midpoint = (x_sum / len(verts), y_sum / len(verts), z_sum / len(verts))
            return midpoint
        midpoints = []
        for face in self.faces:
            curVerts = []
            for f in face:
                curVerts.append(vertices[f])
            midpoints.append(midpoint_of_face(curVerts))

        return midpoints

    def calculate_vertices_and_edges(self) -> list:
        """
        How all the vertices are manually calculated

        :return: List of all vertices defining the 8 corners of the parallelepiped
        """
        v1 = np.array([self.a, 0, 0])
        v2 = np.array([self.b * np.cos(self.gamma), self.b * np.sin(self.gamma), 0])
        v3x = self.c * np.cos(self.beta)
        v3y = self.b * np.cos(self.alpha) * np.sin(self.gamma)
        v3z = np.sqrt(self.c ** 2 - v3x ** 2 - v3y ** 2)
        v3 = np.array([v3x, v3y, v3z])
        vertices = np.around(np.array([
                    np.array([0, 0, 0]),
                    v1,
                    v2,
                    v3,
                    v1 + v3,
                    v2 + v3,
                    v1 + v2,
                    v1 + v2 + v3
                ]), 2)
        return vertices


    def update_face_equations(self) -> None:
        """
        This function will calculate each face equation to use for comparisons and is called when a grammar is applied
        to extend the face directions, for example.
        """
        face_equations = []
        all_vertices = self.calculate_vertices_and_edges()
        # Order of faces matters == LEFT HANDED due to checking for points

        face_vectors = ['xy', 'xy', 'yz', 'yz', 'xz', 'xz']  # Come back to this later, these are the "a b c" plane
        # alignments where x==a y==b and z==c
        for face_vertices in self.faces:
            vertices = []
            for v in face_vertices:
                vertices.append(all_vertices[v])
            vector1 = vertices[1] - vertices[0]
            vector2 = vertices[2] - vertices[0]
            normal = np.cross(vector1, vector2)
            # Ensure the normal points outward from the parallelepiped
            center = np.mean(vertices, axis=0)
            if np.dot(normal, center) > 0:
                normal = -normal
            d = -np.dot(normal, vertices[0])
            a, b, c = normal
            face_equations.append((a, b, c, d))
        self.face_equations = face_equations


    def update_volume(self) -> None:
        """ Calculates and updates the volume attribute """
        self.volume = parallelepiped_volume(self.cell_type, self.a, self.b, self.c, self.alpha,
                                            self.beta, self.gamma)


    @staticmethod
    def intersection_point(line_start: np.array, line_end: np.array, abcd: np.array) -> bool:
        """
        Calculate the intersection point between a line segment and a plane if it exists.

        :param line_start: Starting point of a line segment [X Y Z]
        :param line_end: Ending point of a line segment [X Y Z]
        :param abcd: Defining parameters of a face in Cartesian space

        :return True if there is an intersection, False if not
        """
        # For the current plane we are checking, we use the plane equation and the points of the line. If the value
        # switches from positive to negative (or otherwise) then we must be "crossing" this plane and therefore there
        # is an intersection!
        sign_start = np.dot(abcd[:3], line_start) + abcd[3]
        sign_end = np.dot(abcd[:3], line_end) + abcd[3]

        # We ONLY check in the case that we are not checking from a point that lies on a plane. If a point already
        # lies on a plane then we presume that we are on the correct side.
        if not np.isclose(sign_start, 0) and not np.isclose(sign_end, 0):
            # Check if the line and the plane are parallel and if so just return "None"
            if sign(sign_start) != sign(sign_end):
                return True
            else:
                return False


    def edge_intersecting_box(self, point1: np.ndarray, point2: np.ndarray) -> bool:
        """
        Determines if there is an intersection point between a line segment and a plane in the mesh.

        :param point1: Starting point of a line segment [X Y Z]
        :param point2: Ending point of a line segment [X Y Z]

        :return True if an edge is intersecting a face, false if not
        """
        # Check if a point lies inside the parallelepiped by checking if it lies on the positive side of all faces
        for abcd in self.face_equations:
            # Check if the line segment intersects the face:
            if self.intersection_point(point1, point2, np.array(abcd)):
                return True

        # Otherwise return False meaning "edge not intersecting"
        return False


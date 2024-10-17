"""
This script is specifically used to "hold" the variety of functions used to create a trace for the kodak_creator script
"""
from typing import Optional, Union
import plotly.graph_objs as go
from dataclasses import dataclass
import numpy as np

@dataclass
class PlotTrace(object):
    plotly_trace: Optional[Union[go.Scatter, go.Scatter3d]]  # Note: We update this as we add in new supported data type
    type: str


def scatter(x: Optional[Union[list, np.ndarray]], y: Optional[Union[list, np.ndarray]],
            lines_markers_both: str, **kwargs) -> PlotTrace:
    if lines_markers_both == 'both':
        lines_markers_both = 'markers+lines'
    trace = go.Scatter(x=x, y=y, mode=lines_markers_both)
    for key, value in kwargs.items():
        if hasattr(trace, key):
            trace[key] = value
        else:
            raise Exception(f'Invalid argument for {key} passed in, the go.Scatter object does not support this'
                            f'keyword argument')

    # Now create a PlotTrace to return:
    return PlotTrace(plotly_trace=trace, type='scatter')


def scatter3d(x: Optional[Union[list, np.ndarray]], y: Optional[Union[list, np.ndarray]],
              z: Optional[Union[list, np.ndarray]], **kwargs) -> PlotTrace:

    trace = go.Scatter3d(x=x, y=y, z=z)
    for key, value in kwargs.items():
        if hasattr(trace, key):
            trace[key] = value
        else:
            raise Exception(f'Invalid argument for {key} passed in, the go.Scatter object does not support this'
                            f'keyword argument')

    # Now create a PlotTrace to return:
    return PlotTrace(plotly_trace=trace, type='scatter3d')
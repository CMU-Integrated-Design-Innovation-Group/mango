"""
This script performs a random walk of N steps to validate the grammars and design constraint capabilities. If we were
to get an error, there would be a flaw in the logic somewhere in the code if we consider testing across a large set of
N values.
"""
# Import Modules
from mango.design_spaces.bounding_box import TetragonalBox
from mango.features.preserved_regions import PreservedVertex
from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
from mango.grammars import origami_grammars
from mango.optimization_features import design_constraints
import numpy as np
from copy import deepcopy
import plotly.graph_objs as go

if __name__ == '__main__':
    X = 50
    Y = 50
    Z = 70
    random_seed = 8

    new_box = TetragonalBox(a=X, c=Z)
    preserved_regions = [
        PreservedVertex(v1=np.array([X / 2, 0, Z / 2])),
        PreservedVertex(v1=np.array([X / 2, Y, Z / 2])),
        PreservedVertex(v1=np.array([0, Y / 2, Z / 2])),
        PreservedVertex(v1=np.array([X, Y / 2, Z / 2])),
        PreservedVertex(v1=np.array([X / 2, Y / 2, 0])),
        PreservedVertex(v1=np.array([X / 2, Y / 2, Z])),
    ]

    excluded_regions = []

    design_space = PolyhedralSpace(bounding_box=new_box, preserved=preserved_regions, excluded=excluded_regions)
    design_space.generate_start_shape()

    # Next we import the grammar set and the default design constraints:
    grammar_set = origami_grammars.TriangulationGrammars()
    constraints = design_constraints.PolyhedralDefaultConstraints()
    random_walk_steps = 1000
    rules = grammar_set.grammar_names
    plotting_dict = {}
    cset = ['Outside Design Space', 'Vertex in Excluded', 'Edge in Excluded', 'Invalid Edge Length',
            'Invalid Scaffold Length', 'Invalid Face Angle', 'Broken preserved edge', 'Intersecting edges',
            'Intersecting faces'
            ]
    for r in rules:
        plotting_dict[r] = [0, 0, 0]
    for _ in range(random_walk_steps):
        current_design_state = deepcopy(design_space)
        new_rule = grammar_set.pick_random_grammar()
        grammar_applied_successfully, grammar_time = grammar_set.call_grammar_function(grammar_selected=new_rule,
                                                                                       design_space=current_design_state,
                                                                                       extension_value=0.34)
        constraints.update_params(design_space=current_design_state)
        passed_all_constraints, constraint_times = constraints.check_constraints(design_space=current_design_state,
                                                                                 compare_space=design_space,
                                                                                 constraint_set=cset)
        if grammar_applied_successfully and passed_all_constraints:
            # If the grammar is applied successfully and the constraints are all passed then this new state replaces
            # the current design_space
            del design_space  # Just to be safe (peace of mind) I fully delete design_space from memory and re-assign
            design_space = current_design_state
            plotting_dict[new_rule][0] += 1
        elif grammar_applied_successfully and not passed_all_constraints:
            del current_design_state
            plotting_dict[new_rule][1] += 1
        else:
            plotting_dict[new_rule][2] += 1

    # Plot and show the resultant rule application results:
    groups = list(plotting_dict.keys())
    counts = list(plotting_dict.values())
    name_dict = {
        0: 'Grammar Applied Successfully, passed design constraints',
        1: 'Grammar Applied Successfully, did NOT pass design constraints',
        2: 'Grammar NOT applied Successfully, design constraints irrelevant'
    }
    colors = ['green', 'orange', 'red']
    traces = []
    for i in range(3):  # assuming there are always 3 integers per group
        trace = go.Bar(
            x=groups,
            y=[count[i] for count in counts],
            name=name_dict[i],
            marker=dict(color=colors[i])
        )
        traces.append(trace)
    layout = go.Layout(
        title='Grammar Application Results',
        xaxis=dict(title='Groups'),
        yaxis=dict(title='Count')
    )

    # Create figure
    fig = go.Figure(data=traces, layout=layout)

    # Show figure
    fig.show()


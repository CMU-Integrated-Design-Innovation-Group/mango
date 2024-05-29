# Define the list of specific row indices you care about for each block
import numpy as np
import plotly.graph_objs as go

def process_file(filename, init_node, nodes_to_check):
    current_node_positions = []
    max_distances_per_group = []
    current_data_section = False
    line_counter = -1  # To count the lines within a block
    max_diameters = []

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('t ='):
                # Process the current block before starting a new one:
                if current_node_positions:
                    c1 = check_from_position.split()
                    for i in current_node_positions:
                        c2 = i.split()
                        p1 = np.array([float(c1[0]), float(c1[1]), float(c1[2])])
                        p2 = np.array([float(c2[0]), float(c2[1]), float(c2[2])])
                        max_diameters.append(np.linalg.norm(p2 - p1)*0.8518)  # .8518 is correction factor for oxDNA
                    max_distances_per_group.append(max(max_diameters))
                    max_diameters = []  # Reset

                # Reset for new block
                current_node_positions = []
                check_from_position = None

                line_counter = -1
                current_data_section = False
            elif line.startswith('b =') or line.startswith('E ='):
                current_data_section = True
            else:
                if current_data_section:
                    line_counter += 1
                    if line_counter in nodes_to_check:
                        current_node_positions.append(line)
                    elif line_counter == init_node:
                        check_from_position = line

    # Process final block of time:
    if current_node_positions:
        c1 = check_from_position.split()
        for i in current_node_positions:
            c2 = i.split()
            p1 = np.array([float(c1[0]), float(c1[1]), float(c1[2])])
            p2 = np.array([float(c2[0]), float(c2[1]), float(c2[2])])
            max_diameters.append(np.linalg.norm(p2 - p1)*0.8518)
    max_distances_per_group.append(max(max_diameters))

    return max_distances_per_group

def plot_histo(data, fname):
    plot_layout = go.Layout(
        legend=dict(
            x=0.48,  # Set x position to move the legend to the right
            y=1.1,  # Set y position to move the legend to the top
        ),
        height=300,
        width=300,
        xaxis=dict(
            linecolor='black',
            linewidth=1,
            mirror=False,
            tickfont=dict(size=16, family='arial', color='black'),
            titlefont=dict(size=16, family='arial', color='black'),
            #title='Max hole size (nm)',
            showgrid=False,
            gridcolor='black',
            zeroline=False,
            range=[4, 24]
        ),
        yaxis=dict(
            linecolor='black',
            linewidth=2,
            mirror=False,
            tickfont=dict(size=16, family='arial', color='black'),
            titlefont=dict(size=16, family='arial', color='black'),
            #title='Count',
            showgrid=False,
            gridcolor='black',
            zeroline=False,
            range=[0, 39]
        ),
        font=dict(family='arial', color='black', size=16),
        plot_bgcolor='rgba(0, 0, 0, 0)',
        paper_bgcolor='rgba(0, 0, 0, 0)',
        autosize=False
    )
    fig = go.Figure(data=[go.Histogram(x=data, nbinsx=10, marker=dict(color='orange',  # Bin color
                    line=dict(color='white', width=2)))], layout=plot_layout)
    fig.write_image(fname)


# Specify the filename
design_1 = "./Design1_Output/Design1_trajectory_sim.dat"
design_3 = "./Design3_Output/Design3_trajectory_sim.dat"
design_4 = "./Design4_Output/Design4_trajectory_sim.dat"

# Define out nucleotide #s to track distances of
# These are the indices of the paired nucleotides of each edge at the vertex. We do not consider the unpaired
# nucleotides as they are much more susecptible to noise due to being ssDNA
### initial node needs to be the smallest value (this is so inefficient idk why i coded it like this)
d1_initial_node = 13
d1_hole_nodes_to_check = [1908, 5873, 1914, 5778, 1986, 5773, 1990, 2435, 5756, 5753, 2442, 4163, 2724, 2729, 4159, 5878, 5807,
          5806, 21, 302, 5729, 1737, 5709, 5703, 1741, 4963, 1821, 1828, 4961, 5878]

d3_initial_node = 211
d3_hole_nodes_to_check = [214, 2587, 1072, 2442, 1076, 2436, 1076, 2413, 1173, 2409, 1178, 1240, 3756, 1246, 3752]

d4_initial_node = 393
d4_hole_nodes_to_check = [2218, 2221, 4376, 3472, 4384, 397, 3466, 2025, 4542, 2029, 4536, 4408, 2122, 4402, 2125]


des1_max_hole_size_distribution = process_file(design_1, d1_initial_node, d1_hole_nodes_to_check)

des3_max_hole_size_distribution = process_file(design_3, d3_initial_node, d3_hole_nodes_to_check)

des4_max_hole_size_distribution = process_file(design_4, d4_initial_node, d4_hole_nodes_to_check)

plot_histo(des1_max_hole_size_distribution, 'des1_hole.svg')
plot_histo(des3_max_hole_size_distribution, 'des3_hole.svg')
plot_histo(des4_max_hole_size_distribution, 'des4_hole.svg')





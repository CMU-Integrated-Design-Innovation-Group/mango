# Define the list of specific row indices you care about for each block
import numpy as np
import matplotlib.pyplot as plt
import math

def dist(p1, p2):
    return np.linalg.norm(p2 - p1)

def convert_list(l):
    return_list = []
    for i in l:
        c1 = i.split()
        return_list.append(np.array([float(c1[0]), float(c1[1]), float(c1[2])]))
    return return_list

def find_largest_distance(list1, list2):
    # First convert the lists to a list of points:
    l1 = convert_list(list1)
    l2 = convert_list(list2)

    max_distance = 0

    for point1 in l1:
        for point2 in l2:
            distance = dist(point1, point2)
            if distance > max_distance:
                max_distance = distance

    return max_distance


def process_file(filename, height1, height2, width1, width2, length1, length2):
    current_block_height1 = []
    current_block_height2 = []
    current_block_width1 = []
    current_block_width2 = []
    current_block_length1 = []
    current_block_length2 = []
    current_data_section = False
    line_counter = -1  # To count the lines within a block
    distances = {'height': [],
                 'width': [],
                 'length': []
                 }

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('t ='):
                # Process the current block before starting a new one
                if current_block_height1:
                    for index in ['height', 'width', 'length']:
                        if index == 'height':
                            largest_distance = find_largest_distance(current_block_height1, current_block_height2)
                        elif index == 'width':
                            largest_distance = find_largest_distance(current_block_width1, current_block_width2)
                        else:
                            largest_distance = find_largest_distance(current_block_length1, current_block_length2)

                        distances[index].append(largest_distance*0.8518)  # .8518 is correction factor to nanometers
                # Reset for new block
                current_block_height1 = []
                current_block_height2 = []
                current_block_width1 = []
                current_block_width2 = []
                current_block_length1 = []
                current_block_length2 = []

                line_counter = -1
                current_data_section = False
            elif line.startswith('b =') or line.startswith('E ='):
                current_data_section = True
            else:
                if current_data_section:
                    line_counter += 1
                    if line_counter in height1:
                        current_block_height1.append(line)
                    elif line_counter in height2:
                        current_block_height2.append(line)
                    elif line_counter in width1:
                        current_block_width1.append(line)
                    elif line_counter in width2:
                        current_block_width2.append(line)
                    elif line_counter in length1:
                        current_block_length1.append(line)
                    elif line_counter in length2:
                        current_block_length2.append(line)

    # Process the last block if not empty
    if current_block_height1:
        for index in ['height', 'width', 'length']:
            if index == 'height':
                largest_distance = find_largest_distance(current_block_height1, current_block_height2)
            elif index == 'width':
                largest_distance = find_largest_distance(current_block_width1, current_block_width2)
            else:
                largest_distance = find_largest_distance(current_block_length1, current_block_length2)

            distances[index].append(largest_distance * 0.8518)  # .8518 is correction factor to nanometers

    return distances


def plot_histograms(dists: dict, input_dists: dict):
    plt.rcParams['font.size'] = 18
    plt.rcParams['font.family'] = 'Arial'
    fig, axs = plt.subplots(1, 3, figsize=(14, 4))

    # Plot each histogram
    axs[0].hist(dists['height'], bins=30, color='blue', alpha=0.7)
    axs[0].axvline(input_dists['height'], color='black', linestyle='dashed', linewidth=4)
    axs[0].set_title('Max Height')
    axs[0].set_xlabel('Height (nm)')
    axs[0].set_ylabel('Count')

    axs[1].hist(dists['width'], bins=30, color='green', alpha=0.7)
    axs[1].axvline(input_dists['width'], color='black', linestyle='dashed', linewidth=4)
    axs[1].set_title('Max Width')
    axs[1].set_xlabel('Width (nm)')

    axs[2].hist(dists['length'], bins=30, color='red', alpha=0.7)
    axs[2].axvline(input_dists['length'], color='black', linestyle='dashed', linewidth=4)
    axs[2].set_title('Max Length')
    axs[2].set_xlabel('Length (nm)')

    # Adjust layout
    axs[0].set_xlim([38, 86])
    axs[1].set_xlim([24, 56])
    axs[2].set_xlim([28, 56])
    axs[0].set_ylim([0, 14])
    axs[1].set_ylim([0, 14])
    axs[2].set_ylim([0, 14])
    plt.tight_layout()
    plt.show()


# Specify the filename
design_1 = "./Design1_Output/Design1_trajectory_sim.dat"
design_1_input = '../input_to_oxdna/Design1.dat'

design_3 = "./Design3_Output/Design3_trajectory_sim.dat"
design_3_input = '../input_to_oxdna/Design3.dat'

design_4 = "./Design4_Output/Design4_trajectory_sim.dat"
design_4_input = '../input_to_oxdna/Design4.dat'

# Define out nucleotide #s to track distances of:
d1_h1 = [1343,5180,5179,1349,4869,1138,1139,4862,4840,1193,4833,1196,1297,3164,1301,3162]
d1_h2 = [5228,3036,5203,2826,2829,5196,3773,2880,3769,2886,4251,2969,4247,2973,3030,5230]
d1_w1 = [5073,1686,1684,5079,3269,471,3273,466,3296,411,3299,406,4439,360,4443,357]
d1_w2 = [2384,5467,5471,2380,2335,3858,2331,3863,3885,2277,2274,3890,5620,2042,2041,5625]
d1_l1 = [1828,4960,4963,1821,1908,5878,1914,5873,5702,1742,1737,5709,1986,5778,1990,5773,2435,5756,2442,5753,5729,307,5725,4163,2724,2729,13,5807,5806,21,4159,302]
d1_l2 = [1454,5128,1451,5135,3644,600,605,3640,651,5955,652,5948,1043,3452,1041,3460,6019,949,944,6024,3933,897,3937,895,5301,733,5298,740,826,5906,833,5904]

d3_h1 = [575,4555,2560,692,694,2553,2531,790,2528,794,2847,895,900,4561,570,2844]
d3_h2 = [4119,1365,4125,1361,96,4144,93,4150,3417,13,3420,8,3231,1827,1823,3234]
d3_w1 = [296,3324,391,3303,396,3299,4075,479,483,4071,988,4050,4046,993,292,3328]
d3_w2 = [3902,1452,1456,3898,1551,3877,1735,4189,4194,1732,1637,4215,4219,1633,3872,1554]
d3_l1 = [211,2594,2587,214,1072,2442,2436,1076,2413,1173,2409,1178,1240,3756,3752,1246]
d3_l2 = [1993,4438,1931,3068,4433,1997,4307,2110,2115,4301,4273,2193,2189,4279,3072,1926]

d4_h1 = [21,3862,3843,115,117,3837,210,2603,214,2597,303,4444,4440,308,3868,19]
d4_h2 = [4119,1365,4125,1361,96,4144,93,4150,3417,13,3420,8,3231,1827,1823,3234]
d4_w1 = [1442,2700,1347,2721,1344,2727,4294,1630,1629,4299,1538,3689,1534,3694,2694,1443]
d4_w2 = [2992,678,2473,587,586,2479,491,2500,2986,682,1931,4341,1932,4335,2506,488]
d4_l1 = [3472,393,2221,4376,4383,2218,2125,4402,3466,397,4542,2025,2029,4536,4408,2122]
d4_l2 = [1073,2818,4259,1162,4253,1166,3946,863,866,3939,2837,959,2833,964,1069,2822]


d1_input = process_file(design_1_input, d1_h1, d1_h2, d1_w1, d1_w2, d1_l1, d1_l2)
d3_input = process_file(design_3_input, d3_h1, d3_h2, d3_w1, d3_w2, d3_l1, d3_l2)
d4_input = process_file(design_4_input, d4_h1, d4_h2, d4_w1, d4_w2, d4_l1, d4_l2)
dists1 = process_file(design_1, d1_h1, d1_h2, d1_w1, d1_w2, d1_l1, d1_l2)
dists3 = process_file(design_3, d3_h1, d3_h2, d3_w1, d3_w2, d3_l1, d3_l2)
dists4 = process_file(design_4, d4_h1, d4_h2, d4_w1, d4_w2, d4_l1, d4_l2)

# Create histograms of trajectory:
plot_histograms(dists1, d1_input)
plot_histograms(dists3, d3_input)
plot_histograms(dists4, d4_input)
print(d1_input)
print(d3_input)
print(d4_input)


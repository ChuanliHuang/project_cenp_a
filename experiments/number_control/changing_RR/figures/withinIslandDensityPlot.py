import matplotlib.pyplot as plt
import numpy as np
from itertools import chain
import pandas as pd
import os
import re
import statistics


def show_kymograph(df):
    x = []
    y = []
    for index, row in df.iterrows():
        row = row.to_numpy()
        if row.sum() != 0:
            x_coordinates = np.argwhere(row == 1).flatten().tolist()
            x.append(x_coordinates)
            y_coordinates = [index] * int(row.sum())
            y.append(y_coordinates)
        else:
            break
    x = list(chain.from_iterable(x))
    y = list(chain.from_iterable(y))
    binwidth = 1
    plt.hist2d(y, x,
               bins=[range(min(y), max(y) + binwidth + 1, binwidth), range(min(x), max(x) + binwidth + 1, binwidth)],
               cmap='gray')

    plt.ylabel('nucleosome position')
    plt.xticks(np.arange(0, 2250, step=250), np.arange(0, 112500, step=12500))
    plt.xlabel('generation')
    plt.show()


def get_density_within_islands(df):
    densities = []
    island_nums = []
    for index, row in df.iterrows():
        a = row.to_numpy()
        positions = np.argwhere(a == 1).flatten()
        # deal with periodic boundary
        if positions[0] - (positions[-1] - 5000) < 1000:
            a_ls = np.array_split(a, 2)
            a = np.concatenate((a_ls[1], a_ls[0]), axis=None)
            positions = np.argwhere(a == 1).flatten()
        # recognize islands
        island_threshold = 200
        islands = []
        island = []
        i = 0
        while i < len(positions):
            if i + 1 < len(positions):
                while positions[i + 1] - positions[i] < island_threshold:
                    island.append(positions[i])
                    i += 1
                    if i + 1 >= len(positions):
                        break
            island.append(positions[i])
            islands.append(island)
            island = []
            i += 1
        island_nums.append(len(islands))
        total_island_length = 0
        for j in islands:
            total_island_length += (j[-1] - j[0])
        density = len(positions) / total_island_length
        densities.append(density)
    return densities, island_nums


# loop through sub-folders
# loop through lineages
folder_name = '/Users/kikawaryoku/Desktop/changing_RR/results/'
afs = []
sub_folder_names = []
for sub_folder in os.listdir(folder_name):
    if sub_folder != '.DS_Store':
        sub_folder_name = folder_name + sub_folder + '/states/'
        af = float(sub_folder[-2:]) / 10
        afs.append(af)
        sub_folder_names.append(sub_folder_name)


# zipped = sorted(zip(afs, sub_folder_names))
# for af, sub_folder_name in zipped:
#     print(af)
#     averages = []
#     for lineage in range(10):
#         file_name = sub_folder_name + 'lineage{}.csv'.format(lineage)
#         _df = pd.read_csv(file_name, sep=',', header=None)
#         res, island_num = get_density_within_islands(_df)
#         average = statistics.mean(island_num)
#         averages.append(average)
#     x = af
#     y = statistics.mean(averages)
#     sd = statistics.stdev(averages)
#     plt.errorbar(x, y, yerr=sd, fmt='o', color='black', ecolor='lightgray', elinewidth=3, capsize=0)
# plt.xlabel('af')
# plt.ylabel('average island number')
# plt.show()
#         average = statistics.mean(res)
#         averages.append(average)
#     x = af
#     y = statistics.mean(averages)
#     sd = statistics.stdev(averages)
#     plt.errorbar(x, y, yerr=sd, fmt='o', color='black', ecolor='lightgray', elinewidth=3, capsize=0)
# plt.xlabel('af')
# plt.ylabel('within-island density')
# plt.show()


file_name = '/Users/kikawaryoku/Desktop/changing_RR/results/s3rr3af01/states/lineage1.csv'
_df = pd.read_csv(file_name, sep=',', header=None)
show_kymograph(_df)
res, island_num = get_density_within_islands(_df)
plt.plot(island_num)
plt.ylabel('island number')
plt.ylim(0.5, 8.5)
plt.xlabel('generation')
plt.xticks(np.arange(0, 2250, step=250), np.arange(0, 112500, step=12500))
plt.show()

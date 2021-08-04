import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import re


# loop through sub-folders
# loop through lineages

# af
folder_name = '/Users/kikawaryoku/Desktop/changing_AF/results/'
rrs = []
sub_folder_names = []
pattern = 'rr(.*?)af'
for sub_folder in os.listdir(folder_name):
    if sub_folder != '.DS_Store':
        sub_folder_name = folder_name + sub_folder + '/af_num/'
        substring = re.search(pattern, sub_folder_name).group(1)
        rr = int(substring)
        rrs.append(rr)
        sub_folder_names.append(sub_folder_name)


zipped = sorted(zip(rrs, sub_folder_names))
fig, axs = plt.subplots(4, 1, constrained_layout=True)
j = 0
for rr, sub_folder_name in zipped:
    print(rr)
    frames = []
    for i in range(10):
        file_name = sub_folder_name + 'lineage{}.csv'.format(i)
        _df = pd.read_csv(file_name, sep=',', header=None)
        frames.append(_df)
    # comb_df = pd.concat(frames, axis=1)
    # sd = _df.stack().std()
    # lineage_average = _df.mean()
    # x = af
    # y = lineage_averages
    # plt.errorbar(x, y, yerr=sd, fmt='o', color='black', ecolor='lightgray', elinewidth=3, capsize=0)
    # comb_df_tail = comb_df.tail(1000)
    # rr_num_means = []
    # for index, row in comb_df_tail.iterrows():
    #     rr_num_mean = row.mean()
    #     rr_num_means.append(rr_num_mean)
    if rr == 3 or rr == 4 or rr == 6 or rr == 12:
        axs[j].plot(np.arange(30 * rr) / rr, _df.tail(30 * rr))
        axs[j].plot(np.arange(30 * rr)/rr, _df.tail(30 * rr), '.')
        axs[j].set_title('rr={}'.format(rr))
        axs[j].set_ylabel('af')
        axs[j].set_ylim(-0.05, 1)
        j += 1
plt.xlabel('last 30 generation')
plt.show()

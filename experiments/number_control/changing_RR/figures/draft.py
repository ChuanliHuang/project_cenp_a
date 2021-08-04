import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np


# loop through sub-folders
# loop through lineages

# rr_num
folder_name = '/Users/kikawaryoku/Desktop/changing_RR/results/'
afs = []
sub_folder_names = []
for sub_folder in os.listdir(folder_name):
    if sub_folder != '.DS_Store':
        sub_folder_name = folder_name + sub_folder + '/rr_num/'
        af = float(sub_folder[-2:]) / 10
        afs.append(af)
        sub_folder_names.append(sub_folder_name)


zipped = sorted(zip(afs, sub_folder_names))
for af, sub_folder_name in zipped:
    print(af)
    frames = []
    for i in range(10):
        file_name = sub_folder_name + 'lineage{}.csv'.format(i)
        _df = pd.read_csv(file_name, sep=',', header=None)
        frames.append(_df)
    comb_df = pd.concat(frames, axis=1)
    se = comb_df.stack().sem()
    lineage_average = comb_df.mean().mean()
    x = af
    y = lineage_average
    plt.errorbar(x, y, yerr=se, fmt='o', color='black', ecolor='lightgray', elinewidth=3, capsize=0)
    # comb_df_tail = comb_df.tail(1000)
    # rr_num_means = []
    # for index, row in comb_df_tail.iterrows():
    #     rr_num_mean = row.mean()
    #     rr_num_means.append(rr_num_mean)
#     plt.plot(np.arange(200), _df.tail(200), label='af={}'.format(af))
# plt.legend()
# plt.ylim(-1, 22)
# plt.yticks(np.arange(0, 22, 2))
# plt.xlabel('last 200 generation')
# plt.ylabel('rr number')
plt.show()

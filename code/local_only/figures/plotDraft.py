import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import re

folder_name = '/Users/kikawaryoku/Desktop/local_only/results/'
afs = []
sub_folder_names = []
pattern = 'af(.*?)/'
for sub_folder in os.listdir(folder_name):
    if sub_folder != '.DS_Store':
        sub_folder_name = folder_name + sub_folder + '/states/'
        substring = re.search(pattern, sub_folder_name).group(1)
        while len(substring) < 4:
            substring = substring + '0'
        af = int(substring)
        afs.append(af / 1000)
        sub_folder_names.append(sub_folder_name)

zipped = sorted(zip(afs, sub_folder_names))

x = []
y = []
for af, sub_folder_name in zipped:
    print(af)
    for i in range(10):
        file_name = sub_folder_name + 'lineage{}.csv'.format(i)
        _df = pd.read_csv(file_name, sep=',', header=None)
        x.append(af)
        y.append(_df.mean().mean())

plt.plot(x, y, '.', label='100,000')
plt.xlabel('af')
plt.ylabel('last 1000 density mean')
plt.legend()
plt.show()

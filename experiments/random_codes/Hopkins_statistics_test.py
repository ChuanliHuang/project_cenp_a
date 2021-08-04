from pyclustertend import hopkins
import random
import numpy as np
import pandas as pd

'''A value for H higher than 0.75 indicates a clustering tendency at the 90% confidence level.'''


def init_gen_0(nuc_num, init_den):
    gen_0_nuc_comp = np.zeros(nuc_num, dtype=int)
    for i in range(nuc_num):
        if random.random() < init_den:
            gen_0_nuc_comp[i] = 1
    return gen_0_nuc_comp


random_array = init_gen_0(1000, 0.2)
x = np.argwhere(random_array == 1).flatten()
print(hopkins(x, len(x)))

file_name = '/experiments/local_only/master/states/lineage0.csv'
df = pd.read_csv(file_name, sep=',', header=None, dtype='int64')
clustered_array = df.loc[20, :].values.flatten()
x = np.argwhere(clustered_array == 1).flatten()
print(hopkins(x, len(x)))
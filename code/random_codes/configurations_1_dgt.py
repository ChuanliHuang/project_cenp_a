import numpy as np
# configurations
trial = 2
nucleosome_num = 2501
cenp_a_fraction = 0.04
dist_bias = 0.5
generations_to_simulate = 10000
mu, sigma = 0, 1.4  # mean and standard deviation
replenish_rounds = 3
amplitude_factor = 0.7
amplitude_factor_long = 0.015  # 0.015
initiate_cell_number = 10
dt_1_dgt = np.dtype([('generation', np.int32), ('lineage', np.int32),
               ('nucleosome_composition', np.int32, (nucleosome_num,))])  # defining data type of structured arrays
data = np.array([], dtype=dt_1_dgt)  # an array to store structured arrays
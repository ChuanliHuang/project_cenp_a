import numpy as np
# configurations
nucleosome_num = 2500
cenp_a_fraction = 0.04
dist_bias = 0.5
generations_to_simulate = 13
replenish_rounds = 3
amplitude_factor = 1
initiate_cell_number = 1
density_max_variation = 100
mu, sigma = 0, 1
dt = np.dtype([('cell_id', np.int32), ('generation', np.int32), ('parent', np.int32),
               ('nucleosome_composition', np.int32, (nucleosome_num,))])  # defining data type of structured arrays
data = np.array([], dtype=dt)  # an array to store structured arrays

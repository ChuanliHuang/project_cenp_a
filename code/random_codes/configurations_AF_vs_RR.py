import numpy as np
# configurations
nucleosome_num = 2500
cenp_a_fraction = 0.04
dist_bias = 0.5
generations_to_simulate = 1000
mu= 0
sigma = 3
initiate_cell_number = 1
af_start = 0
af_stop = 1
af_step = 101
RR_start = 2
RR_stop = 10
RR_step = 9
dt_1_dgt = np.dtype([('generation', np.int32), ('lineage', np.int32),
               ('nucleosome_composition', np.int32, (nucleosome_num,))])
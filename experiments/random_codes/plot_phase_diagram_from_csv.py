import numpy as np
import matplotlib.pyplot as plt
from configurations_1_dgt_new import *
file_name = '/Users/kikawaryoku/PycharmProjects/project_cenp_a/exp_3/raw_data/data_1000.csv'
arr = np.loadtxt(file_name, delimiter=',')
plt.imshow(arr, cmap='cividis')
plt.colorbar()
x_ls = np.linspace(sigma_start, sigma_stop, 11).round(2).tolist()
y_ls = np.linspace(af_start, af_stop, 11).round(2).tolist()
x_ticks = map(str, x_ls)
y_ticks = map(str, y_ls)
plt.xlabel('sigma')
plt.xticks(np.linspace(0, sigma_step - 1, 11), x_ticks)
plt.yticks(np.linspace(0, af_step - 1, 11), y_ticks)
plt.ylabel('loading efficiency')
plt.title('1000-1200')
plt.show()
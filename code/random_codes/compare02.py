import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


file_name = '/0_code/stochastic_CA/master/states/lineage0.csv'
df = pd.read_csv(file_name, sep=',', header=None, dtype='int64')
x = []
y = []
row = df.iloc[0]
row = row.to_numpy()
print(row.sum() / 5000)
x_coordinates = np.argwhere(row == 1).flatten().tolist()
y = [1] * int(row.sum())
plt.plot(x_coordinates, y, '.', markersize=0.5, label='ρ={}'.format(str(row.sum() / 5000)))



file_name = '/0_code/stochastic_CA/master/start_from_005/lineage0.csv'
df = pd.read_csv(file_name)
x = []
y = []
row = df.iloc[2]
row = row.to_numpy()
x_coordinates = np.argwhere(row == 1).flatten().tolist()
y = [0] * int(row.sum())
plt.plot(x_coordinates, y, '.', markersize=0.5, label='ρ={}'.format(str(row.sum() / 5000)))


plt.yticks([0, 1], ['ρ$_{0}$=0.05', 'ρ$_{0}$=0.2'])
plt.legend()
plt.ylim(-1, 2)
plt.xlabel('positions')
plt.show()
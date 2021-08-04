import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

file_name = '/Users/kikawaryoku/OneDrive - University of Edinburgh/CENP-A/important_data/3_3_AF_test.xlsx'
df = pd.read_excel(file_name, engine='openpyxl', sheet_name=2)
plt.plot(df.AF, df.density, '.', label='1')



file_name = '/0_code/stochastic_CA/analysis/bifurcation/means/last_1000_means_al=5000_old.csv'
df = pd.read_csv(file_name, index_col=0)
alphas = []
densities = []
for i, row in df.iterrows():
    alphas.append(i)
    densities.append(row.mean())
plt.plot(alphas, densities, '.', label='2')


plt.legend()
plt.xlabel('af')
plt.ylabel('d_100,000')
plt.show()

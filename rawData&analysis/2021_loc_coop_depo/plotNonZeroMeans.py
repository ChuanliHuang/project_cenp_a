import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import sem

file_name = 'init_den=005/means.csv'
df = pd.read_csv(file_name, index_col=0)
for index, row in df.iterrows():
    row = row[row != 0]
    if len(row) > 0:
        a = plt.errorbar(index, row.mean(), marker='o', yerr=sem(row), color='tab:blue')

file_name = 'init_den=05/means.csv'
df = pd.read_csv(file_name, index_col=0)
for index, row in df.iterrows():
    row = row[row != 0]
    row = row[row != 0]
    if len(row) > 0:
        b = plt.errorbar(index, row.mean(), marker='o', yerr=sem(row), color='tab:orange')
plt.legend([a, b], [r'$\rho_{initial}=0.05$', r'$\rho_{initial}=0.5$'])
plt.title(r'non-zero means')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'asymptotic $\rho$')
plt.show()

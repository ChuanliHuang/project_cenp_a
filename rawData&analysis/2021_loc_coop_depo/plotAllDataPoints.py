import pandas as pd
import matplotlib.pyplot as plt

file_name = 'init_den=005/means.csv'
df = pd.read_csv(file_name, index_col=0)
for index, row in df.iterrows():
    for lineage in row:
        a, = plt.plot(index, lineage, 'o', color='tab:blue')
file_name = 'init_den=05/means.csv'
df = pd.read_csv(file_name, index_col=0)
for index, row in df.iterrows():
    for lineage in row:
        b, = plt.plot(index, lineage, '.', color='tab:orange')
plt.legend([a, b], [r'$\rho_{initial}=0.05$', r'$\rho_{initial}=0.5$'])
plt.title(r'$\sigma =3 \ \omega =7$')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'asymptotic $\rho$')
plt.show()

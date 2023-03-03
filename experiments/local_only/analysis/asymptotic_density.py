import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['savefig.dpi'] = 300

file_name = '/Users/kikawaryoku/Library/CloudStorage/OneDrive-UniversityofEdinburgh/CENP-A/important_data/3_3_AF_test.xlsx'
df_1 = pd.read_excel(file_name, sheet_name=0)
df_2 = pd.read_excel(file_name, sheet_name=2)
plt.plot(df_1.AF, df_1.density, '.', label='10,000 step simulation')
plt.plot(df_2.AF, df_2.density, '.', label='100,000 step simulation')
plt.ylabel(r'$O(\rho$)')
plt.xlabel(r'$\alpha$')
plt.legend()
plt.show()

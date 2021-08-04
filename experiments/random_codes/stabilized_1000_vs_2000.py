import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
file_name = '/Users/kikawaryoku/PycharmProjects/project_cenp_a/exp_4/raw_data/big_table.csv'
df = pd.read_csv(file_name)
# for i in df['S'].unique():
#     df_S = df[df['S'] == i]
#     x = df_S['AF'].unique()
#     y = []
#     z = []
#     for j in df_S['AF'].unique():
#         df_S_AF = df_S[df_S['AF'] == j]
#         try:
#             true_count_1000 = df_S_AF['stabilized_1000'].value_counts()[True]
#         except:
#             true_count_1000 = 0
#         try:
#             true_count_2000 = df_S_AF['stabilized_2000'].value_counts()[True]
#         except:
#             true_count_2000 = 0
#         y.append(true_count_1000)
#         z.append(true_count_2000)
#     plt.title('RR={} S={}'.format(3, i))
#     plt.ylabel('number of \'stabilized\' lineages')
#     plt.xlabel('loading efficiency')
#     plt.plot(x, y, label='1000')
#     plt.plot(x, z, label='2000')
#     plt.legend()
#     plt.show()

for i in df['S'].unique():
    df_S = df[df['S'] == i]
    x = df_S['AF'].unique()
    y = []
    z = []
    y_err = []
    z_err = []
    for j in df_S['AF'].unique():
        df_S_AF = df_S[df_S['AF'] == j]
        selection = df_S_AF['stabilized_2000']
        selected_df = df_S_AF[selection]
        selected_df = df_S_AF
        std_1000 = selected_df['mean_1000'].std()
        std_2000 = selected_df['mean_2000'].std()
        if len(selected_df) > 0:
            mean_1000 = selected_df['mean_1000'].mean()
            mean_2000 = selected_df['mean_2000'].mean()
        else:
            mean_1000 = 0
            mean_2000 = 0
        y.append(mean_1000)
        y_err.append(std_1000)
        z.append(mean_2000)
        z_err.append(std_2000)
    plt.errorbar(x, y, yerr=y_err, capsize=4, label='1000')
    plt.errorbar(x, z, yerr=z_err, capsize=4, label='2000')
    plt.title('RR={} S={}'.format(3, i))
    plt.ylabel('mean CENP-A density')
    plt.xlabel('loading efficiency')
    plt.legend()
    plt.show()
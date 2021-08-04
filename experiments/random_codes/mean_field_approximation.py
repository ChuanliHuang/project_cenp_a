import math
import numpy as np
import matplotlib.pyplot as plt


def non_zero_probability(_s):
    return 1 - math.erf(1 / (_s * 2 ** 1.5))


def f(_x):
    return _x + _x * af * non_zero_probability(s) * (1 - _x)

def g(_x):
    for i in range(0, rr):
        _x = f(_x)
    return (_x + noise_plus * (1 - _x) - noise * _x * (1 - _x)) * 0.5


rr = 3
af = 1
s = 3
noise_plus = 0
noise = 0.008

def initialize():
    global x, res
    res = []
    x = 0.02

def observe():
    global x, res
    res.append(x)

def update():
    global x, res
    x = g(x)

# initialize()
# for t in range(30):
#     observe()
#     update()


x = np.arange(0, 1.001, 0.001)
y = [i for i in x]
y = np.array(y)
plt.plot(x, y, '--', color='k')
for i in np.arange(0, 1.2, 0.2):
    noise = i
    z = [g(i) for i in x]
    z = np.array(z)
    idx = np.argwhere(np.diff(np.sign(y - z))).flatten()
    plt.plot(x, z, label='noise={}'.format(round(noise, 2)))
    plt.plot(x[idx], y[idx], 'r.')
plt.legend()
plt.xlabel('d_t')
plt.ylabel('d_t+1')
plt.title('af={} s={} rr={}'.format(af, s, rr))
plt.show()

# x = np.arange(0, 1, 0.001)
# y = [i for i in x]
# y = np.array(y)
# prediction = []
# for i in np.arange(0, 1, 0.02):
#     af = i
#     z = [g(i) for i in x]
#     z = np.array(z)
#     idx = np.argwhere(np.diff(np.sign(y - z))).flatten()
#     prediction.append(x[idx][-1])
# plt.plot(np.arange(0, 1, 0.02), prediction, '.', color='royalblue', label='mean-field approximation', zorder=10)
#
# import pandas as pd
#
# file_name = '/Users/kikawaryoku/Desktop/3_3_AF_test.xlsx'
# df = pd.read_excel(file_name, engine='openpyxl', sheet_name=0)
# x = np.concatenate((np.arange(0, 0.4, 0.01), df.AF.to_numpy()))
# y = np.concatenate((np.zeros(40), df.density.to_numpy()))
# plt.plot(x, y, '.', color='tab:red', label='10,000-step simulation', zorder=5)
#
# df = pd.read_excel(file_name, engine='openpyxl', sheet_name=1)
# plt.plot(df.AF, df.density, '.', color='tab:orange', label='50,000-step simulation')
#
# df = pd.read_excel(file_name, engine='openpyxl', sheet_name=2)
# plt.plot(df.AF, df.density, '.', color='tab:green', label='100,000-step simulation')
#
# plt.title('s=3 rr=3')
# plt.xlabel('af')
# plt.ylabel('d_final')
# plt.legend()
# plt.show()





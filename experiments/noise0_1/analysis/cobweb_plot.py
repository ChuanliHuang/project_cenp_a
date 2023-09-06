import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

alpha = 1
s = 3
rr = 1
noise = 0

def non_zero_probability(_s):
    return 1 - math.erf(1 / (_s * 2 ** 1.5))


def f(_x):
    return _x + (1 - _x) * non_zero_probability(s) * alpha * _x


def g(_x):
    for i in range(0, rr):
        _x = f(_x)
    return (_x + noise * (1 - _x)) * 0.5


def func(_x):
    return [_x[1] - g(_x[0]), _x[1] - _x[0]]


def initialize(d_0):
    global data, density
    data = []
    density = d_0
    data.append(density)


def update():
    global data, density
    density = f(density)


def observe():
    global data, density
    data.append(density)

# approaching steady state
# for i in np.arange(0.05, 0.55, 0.05):
#     initialize(i)
#     for j in range(50):
#         update()
#         observe()
#     data = [x - 0.42 for x in data]
#     data = np.array(data)
#     data = np.absolute(data)
#     data = np.log(data)
#     plt.plot(data, label='d0={}'.format(str(round(i, 2))))
# plt.ylim(-12, 0)
# plt.xlim(-0.5, 50.5)
# plt.ylabel('ln |d_0 - d_steady|')
# plt.xlabel('generation')
# plt.title(r'$\alpha={}$'.format(alpha))
# plt.legend()
# plt.show()


x = np.arange(0, 1.0001, 0.0001)
y = [i for i in x]
y = np.array(y)
plt.plot(x, y, '--', color='k', linewidth=0.5)
# for i in np.arange(0, 0.12, 0.02):
for i in [0, 0.001, 0.01, 0.1, 1]:
    noise = i
    z = [g(i) for i in x]
    z = np.array(z)
    plt.plot(x, z, linewidth=0.5, label=r'$noise_{0\rightarrow1}$' + '={}'.format(round(noise, 3)))
    root = optimize.fsolve(func, [0.5, 0.5])
    plt.plot(root[0], root[1], '.', markersize=3, color='black', fillstyle='none')
# plt.legend()
plt.xlim(-0.005, 0.105)
plt.ylim(-0.005, 0.105)
plt.xlabel(r'$\rho_{t}$')
plt.ylabel(r'$\rho_{t+1}$')
# plt.title(r'$\alpha={}$ $rr={}$'.format(alpha, rr))
plt.show()


# x = np.arange(0, 1, 0.001)
# y = [i for i in x]
# y = np.array(y)
# plt.plot(x, y, '--', color='k')
# for i in np.arange(0, 1.2, 0.2):
#     beta = i
#     z = [f(i) for i in x]
#     z = np.array(z)
#     idx = np.argwhere(np.diff(np.sign(y - z))).flatten()
#     plt.plot(x, z, label='β={}'.format(round(beta, 2)))
#     plt.plot(x[idx], y[idx], 'r.')
# plt.legend()
# plt.xlabel('d_t')
# plt.ylabel('d_t+1')
# plt.title('α={}'.format(alpha))
# plt.show()

# x = np.arange(0, 1, 0.001)
# y = [i for i in x]
# y = np.array(y)
# prediction = []
# for i in np.arange(0, 5, 0.05):
#     alpha = i
#     z = [f(i) for i in x]
#     z = np.array(z)
#     idx = np.argwhere(np.diff(np.sign(y - z))).flatten()
#     prediction.append(x[idx][-1])
# plt.plot(np.arange(0, 5, 0.05), prediction, '.', color='royalblue', label='mean-field approximation', zorder=10)
# plt.ylim(-0.05, 0.55)
# plt.title('β={}'.format(beta))
# plt.xlabel('α')
# plt.ylabel('D_stable')
# plt.show()
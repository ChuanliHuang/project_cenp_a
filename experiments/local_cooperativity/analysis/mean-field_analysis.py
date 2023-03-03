import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


def non_zero_probability(_s):
    return 1 - math.erf(1 / (_s * 2 ** 1.5))


def f(_x):
    # return _x + _x * (1 - (1-_x)**2) * af * non_zero_probability(s) * (1 - _x)
    return _x + _x**2 * af * non_zero_probability(s) * (1 - _x)


def g(_x):
    for i in range(0, rr):
        _x = f(_x)
    return _x * 0.5


def func(_x):
    return [_x[1] - g(_x[0]), _x[1] - _x[0]]


rr = 5
af = 1
s = 3


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


x = np.arange(0, 1.001, 0.001)
y = [i for i in x]
y = np.array(y)
plt.plot(x, y, '--', color='k', linewidth=0.5)
for i in np.arange(0, 1.2, 0.2):
    af = i
    z = [g(i) for i in x]
    z = np.array(z)
    plt.plot(x, z, linewidth=0.5, label=r'$\alpha$={}'.format((round(af, 2))))
    root = optimize.fsolve(func, [0.25, 0.25])
    if root[0] == root[1]:
        plt.plot(root[0], root[1], '.', markersize=3, color='black', fillstyle='none')
    root = optimize.fsolve(func, [0.5, 0.5])
    if root[0] == root[1]:
        plt.plot(root[0], root[1], '.', markersize=3, color='black', fillstyle='none')
plt.legend()
plt.xlabel(r'$\rho_{t}$')
plt.ylabel(r'$\rho_{t+1}$')
# plt.title(r'$\sigma={}$ $\omega={}$'.format(s, rr))
plt.show()

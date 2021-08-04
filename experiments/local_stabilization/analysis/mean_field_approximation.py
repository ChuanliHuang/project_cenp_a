import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


def non_zero_probability(_s):
    return 1 - math.erf(1 / (_s * 2 ** 1.5))


def f(_x):
    return _x + _x * alpha * non_zero_probability(s) * (1 - _x)


def h(_x):
    return _x - _x * beta * (1 - _x)  # cooperating stabilization
    # return _x - _x * beta * non_zero_probability(s) * (1 - _x)  # local stabilization
# def h(_x_old, _x):
#     return _x - beta * (_x - _x_old) * (1 - _x_old)  # only new needs stabilization


def g(_x):
    _x_old = _x
    for i in range(0, rr):
        _x = f(_x)
    for j in range(0, 1):
        _x = h(_x)
    return (_x + noise) * 0.5


def func(_x):
    return _x - g(_x)


rr = 3
alpha = 1
s = 3
beta = 1
noise = 0


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


x = np.arange(0, 1.0001, 0.0001)
y = [i for i in x]
y = np.array(y)
plt.plot(x, y, '--', color='k', linewidth=0.5)
# for i in np.arange(0, 0.12, 0.02):
for i in [0, 0.2, 0.4, 0.6, 0.8, 1]:
    beta = i
    z = [g(i) for i in x]
    z = np.array(z)
    plt.plot(x, z, linewidth=0.5, label=r'$\beta$'+'={}'.format(round(beta, 2)))
    # root = optimize.fsolve(func, [0, 0.005, 0.1])
    # plt.plot(root, root, '.', markersize=3, color='black', fillstyle='none')
plt.legend()
plt.xlabel(r'$\rho_{t}$')
plt.ylabel(r'$\rho_{t+1}$')
plt.title(r'$\alpha$={} $s={}$ $rr={}$'.format(alpha, s, rr))
plt.show()





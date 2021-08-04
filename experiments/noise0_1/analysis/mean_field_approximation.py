import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


def non_zero_probability(_s):
    return 1 - math.erf(1 / (_s * 2 ** 1.5))


def f(_x):
    return _x + _x * af * non_zero_probability(s) * (1 - _x)


def g(_x):
    for i in range(0, rr):
        _x = f(_x)
    return (_x + noise * (1 - _x)) * 0.5


def func(_x):
    return [_x[1] - g(_x[0]), _x[1] - _x[0]]


rr = 1
af = 1
s = 3
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
for i in [0, 0.01, 0.1, 1]:
    noise = i
    z = [g(i) for i in x]
    z = np.array(z)
    plt.plot(x, z, linewidth=0.5, label=r'$noise_{0\rightarrow1}$'+'={}'.format(round(noise, 2)))
    root = optimize.fsolve(func, [0, 0])
    plt.plot(root[0], root[1], '.', markersize=3, color='black', fillstyle='none')
plt.legend()
plt.xlabel(r'$\rho_{t}$')
plt.ylabel(r'$\rho_{t+1}$')
plt.title(r'$af={}$ $s={}$ $rr={}$'.format(af, s, rr))
plt.show()




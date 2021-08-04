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
    if _x > threshold * 2:
        return (_x + noise * (1 - _x)) * 0.5
    else:
        return noise * (1 - _x) * 0.5

def func(_x):
    return [_x[1] - g(_x[0]), _x[1] - _x[0]]


rr = 3
af = 0.295
s = 3
noise = 0 # 0.001
threshold = 0# 0.01

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
for i in [0, 0.2, 0.4, 0.6, 0.8, 1]:
    af = i
    z = [g(i) for i in x]
    z = np.array(z)
    plt.plot(x, z, linewidth=0.5, label=r'$\alpha$'+'={}'.format(af))
    root = optimize.fsolve(func, [0.5, 0.5])
    plt.plot(root[0], root[1], '.', markersize=3, color='black', fillstyle='none')
# z = [g(i) for i in x]
# z = np.array(z)
# plt.plot(x, z, linewidth=0.5)
# root = optimize.fsolve(func, [0.05, 0.05])
# plt.plot(root, root, '.', markersize=3, color='black', fillstyle='none')
# root = optimize.fsolve(func, [0, 0])
# plt.plot(root, root, '.', markersize=3, color='black', fillstyle='none')
# plt.xlim(-0.005, 0.045)
# plt.ylim(-0.005, 0.045)
plt.xlabel(r'$\rho_{t}$')
plt.ylabel(r'$\rho_{t+1}$')
plt.legend()
# plt.title(r'$\alpha={}$ $s={}$ $rr={}$ $noise0\rightarrow1={}$ $threshold={}$'.format(af, s, rr, noise, threshold))
plt.show()





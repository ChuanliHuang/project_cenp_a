import matplotlib.pyplot as plt
import numpy as np
import math

def initiate(x0):
    global x
    global res
    x = x0
    res = [x]

def observe():
    global res
    res.append(x)

def update():
    global x
    # x = 0.35 + x
    # n = 4
    # x = 0.35 * x**n / (0.175**n + x**n) + x
    x = 0.6 * (1 - math.e ** (-2.5 * x)) + x
    # x = A * math.log(B * x + 1, 10) + x
    if x > 1:
        x = 1
    x = 0.5 * x


for i in np.linspace(0.05, 1, 20):
    initiate(i)
    for j in range(30):
        update()
        observe()
    res = np.array(res)
    res = np.around(res, 5)
    plt.plot(res)
plt.title('discrete-time model')
plt.xlabel('t')
# plt.xlim(-1, 5)
plt.ylabel('x')
# plt.semilogy()
plt.show()
import numpy as np
import matplotlib.pyplot as plt

eta = 0.2

def f(m):
    return -2*m**3+2*(1+eta)*m**2-2*(1/F+eta)*m+1/F
Fs = np.linspace(1, 20, 100)
for F in Fs:
    x = np.linspace(0, 1, 1000)
    y = [f(i) for i in x]
    z = [0 for i in x]
    y = np.array(y)
    z = np.array(z)
    idx = np.argwhere(np.diff(np.sign(y - z))).flatten()
    plt.plot([F for i in x[idx]], x[idx], '.', color='black', markersize=3)
plt.xlim(1, 20)
plt.ylim(0, 1)
plt.title(r'$\eta$={}'.format(eta))
plt.xlabel('$F$')
plt.ylabel('$m_{eq}$')
plt.show()

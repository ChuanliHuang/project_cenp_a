import math
import numpy as np
import matplotlib.pyplot as plt


N = 5000
x = np.arange(0, N+1, 1)
y = [math.log10(math.comb(N, i)) for i in x]
x = x / N
plt.plot(x, y)
plt.xlabel('ρ')
plt.ylabel('log' + r'$_{10}$' + '$\Omega$')
plt.title('N={}'.format(N))
plt.show()
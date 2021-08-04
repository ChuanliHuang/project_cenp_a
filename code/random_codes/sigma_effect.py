import numpy as np
import matplotlib.pyplot as plt
for sigma in range(3, 11):
    mu = 0
    size = 10000
    rn = np.random.normal(mu, sigma, size=size)
    rn = np.round(rn)
    bars, height = np.unique(rn, return_counts=True)
    density = height / size
    plt.bar(bars, density)
    plt.title('sigma = {}'.format(round(sigma, 1)))
    plt.ylabel('density')
    plt.xlabel('position')
    plt.show()
    plt.clf()
import matplotlib.pyplot as plt
import math
import numpy as np
x = np.linspace(0, 1)
A = 0.52
B = 10
y = []
for i in x:
    y_i = A * math.log(B * i + 1, 10) + i
    if y_i > 1:
        y_i = 1
    y_i = 0.5 * y_i
    y.append(y_i)
plt.plot(x, y)
plt.plot(x, x, '--')
plt.show()

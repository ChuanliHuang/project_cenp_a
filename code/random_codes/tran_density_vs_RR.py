import matplotlib.pyplot as plt
x = range(2, 11)
y = [0.154, 0.152, 0.182, 0.187, 0.2, 0.194, 0.173, 0.193, 0.157]
plt.scatter(x, y)
plt.xlabel('replenishment rounds')
plt.ylabel('CENP-A density at transition')
plt.title('S=3')
plt.show()
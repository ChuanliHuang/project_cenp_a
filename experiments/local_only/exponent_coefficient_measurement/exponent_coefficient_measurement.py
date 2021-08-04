import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model
from sklearn.metrics import r2_score
import math

# file_name = '/Users/kikawaryoku/Desktop/3_3_AF_test.xlsx'
# df = pd.read_excel(file_name, engine='openpyxl', sheet_name=0)
# selection = df.AF > 0.56
# selected_df = df[selection]
# af = selected_df.AF - 0.53
# af = np.log(af)
# af = af.to_numpy()
# af = af.reshape(-1, 1)
# density = selected_df.density
# density = np.log(density)


def non_zero_probability(_s):
    return 1 - math.erf(1 / (_s * 2 ** 1.5))


def f(_x):
    return _x + _x * af * non_zero_probability(s) * (1 - _x)

def g(_x):
    for i in range(0, rr):
        _x = f(_x)
    return _x * 0.5


rr = 3
af = 1
s = 3
x = np.arange(0, 1, 0.001)
y = [i for i in x]
y = np.array(y)
prediction = []
for j in np.arange(0.33, 1, 0.01):
    af = j
    z = [g(i) for i in x]
    z = np.array(z)
    idx = np.argwhere(np.diff(np.sign(y - z))).flatten()
    prediction.append(x[idx][-1])

af = np.arange(0.33, 1, 0.01)
af -= 0.3
af = np.log(af)
af = af.reshape(-1, 1)
density = np.array(prediction)
density = np.log(density)


# Create linear regression object
regr = linear_model.LinearRegression()

# Train the model
regr.fit(af, density)

# Predicition
y_pred = regr.predict(af)

slope = regr.coef_
r2 = r2_score(density, y_pred)
a = np.linspace(-5, 0, 100)
b = slope[0] * a + regr.intercept_

plt.plot(af, density, '.', color='royalblue', label='mean-field approximation')
plt.plot(a, b, '--', color='grey', label='fit')
plt.xlabel('log(af - af_c)')
plt.ylabel('log(d_final)')
slope = round(slope[0], 2)
r2 = round(r2, 3)
plt.title('a={} r^2={}'.format(slope, r2))
plt.legend()
plt.show()

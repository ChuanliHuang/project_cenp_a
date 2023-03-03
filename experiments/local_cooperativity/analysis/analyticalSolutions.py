from sympy import Symbol
from sympy import Interval
from sympy import solveset, S
import math


def non_zero_probability(_s):
    return 1 - math.erf(1 / (_s * 2 ** 1.5))


def f(_x):
    # return _x + _x * (1 - (1-_x)**2) * af * non_zero_probability(s) * (1 - _x)
    return _x + _x ** 2 * af * non_zero_probability(s) * (1 - _x)


def g(_x):
    for i in range(rr):
        _x = f(_x)
    return _x


s = 3
rr = 3
af = 1

af = Symbol('af')
x = Symbol('x')
sol = solveset(g(x) / 2 - x, x, domain=Interval(0, 0.5))
print(sol)

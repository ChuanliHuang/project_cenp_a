import random
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import pycxsimulator

nucleosomeNumber = 5000
initialDensity = 0.02
dilutionBias = 0.5
sigma = 1  # standard deviation of weight function (Gaussian)
alpha = 3.5  # arbitrary parameter
beta = 0  # loss rate


def initial_density(val=initialDensity):
    global initialDensity
    initialDensity = float(val)
    return val


def number_of_nucleosomes(val=nucleosomeNumber):
    global nucleosomeNumber
    nucleosomeNumber = int(val)
    return val


def dilution_bias(val=dilutionBias):
    global dilutionBias
    dilutionBias = float(val)
    return val


def sigma_val(val=sigma):
    global sigma
    sigma = int(val)
    return val


def alpha_val(val=alpha):
    global alpha
    alpha = float(val)
    return val


def beta_val(val=beta):
    global beta
    beta = float(val)
    return val


def generate_discrete_gaussian_distribution(s):
    """ Pre-calculate the probabilities based on normal distribution
      return the generating function
      Wings that further than 3*sigma will be rounded
  """
    left_wing = int(np.round(0 - 3 * s))
    right_wing = int(np.round(0 + 3 * s))
    shifts = np.arange(left_wing, right_wing + 1)
    # get probability for each integer within shifts
    probs = stats.norm.cdf(shifts + 0.5, 0, s) - stats.norm.cdf(shifts - 0.5, 0, s)
    # we need to add all the wings to the outermost bins
    probs[0] = stats.norm.cdf(left_wing + 0.5, 0, s)
    probs[-1] = 1 - stats.norm.cdf(right_wing - 0.5, 0, s)
    # check that sum of probabilities is 1
    assert (np.isclose(sum(probs), 1))
    return probs


def init_gen_0(nuc_num, init_den):
    gen_0_state = np.zeros(nuc_num, dtype=int)
    for i in range(nuc_num):
        if random.random() < init_den:
            gen_0_state[i] = 1
    return gen_0_state


def replenish(old_state, s, a, b):
    new_state = np.copy(old_state)
    zero_positions = np.argwhere(old_state == 0).flatten()
    one_positions = np.argwhere(old_state == 1).flatten()
    f = generate_discrete_gaussian_distribution(s)
    # incorporation of new CENP-A
    for x in zero_positions:
        p_incorp = 0
        for dx in range(3 * -s, 3 * s + 1):
            p_incorp += old_state[((x + dx) % len(old_state))] * f[dx + 3 * s]
        p_incorp *= a
        print(p_incorp)
        if random.random() < p_incorp:
            new_state[x] = 1
    # loss of old CENP-A
    for y in one_positions:
        if random.random() < b:
            new_state[y] = 0
    return new_state


def dilute(old_state, bias):
    new_state = np.copy(old_state)
    one_positions = np.argwhere(old_state == 1).flatten()
    for i in one_positions:
        if random.random() < bias:
            new_state[i] = 0
    return new_state


def initialize():
    global state
    state = init_gen_0(nucleosomeNumber, initialDensity)


def observe():
    global state
    density = state.sum() / len(state)
    positions = np.argwhere(state == 1).flatten()
    plt.cla()
    theta = 2 * np.pi * positions / len(state)
    r = [1] * len(positions)
    plt.polar(theta, r, '.')
    plt.title('density={}'.format(density))
    plt.axis('off')


def update():
    global state
    state = replenish(state, sigma, alpha, beta)
    state = dilute(state, dilutionBias)


pycxsimulator.GUI(
    parameterSetters=[initial_density, number_of_nucleosomes, dilution_bias, sigma_val, alpha_val, beta_val]).start(
    func=[initialize, observe, update])

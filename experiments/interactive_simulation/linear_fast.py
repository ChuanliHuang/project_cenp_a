import scipy.stats as stats
import numpy as np
import random
import matplotlib.pyplot as plt
import pycxsimulator

nucleosome_num = 5000
cenp_a_fraction = 0.02
dist_bias = 0.5
mu = 0
sigma = 3
replenish_rounds = 3
amplitude_factor = 0.5

def loading_efficiency(val=amplitude_factor):
    global amplitude_factor
    amplitude_factor = float(val)
    return val


def rounds_of_replenishment(val=replenish_rounds):
    global replenish_rounds
    replenish_rounds = int(val)
    return val


def loading_range(val=sigma):
    global sigma
    sigma = float(val)
    return val


def initial_density(val=cenp_a_fraction):
    global cenp_a_fraction
    cenp_a_fraction = float(val)
    return val


def number_of_nucleosomes(val=nucleosome_num):
    global nucleosome_num
    nucleosome_num = int(val)
    return val


def distribution_bias(val=dist_bias):
    global dist_bias
    dist_bias = float(val)
    return val


def generate_probability_function(mu, sigma, amplitude_factor):
    """ Pre-calculate the probabilities based on normal distribution
      return the generating function
      Wings that further than 3*sigma will be rounded
  """
    left_wing = int(np.round(mu - 3 * sigma))
    right_wing = int(np.round(mu + 3 * sigma))
    shifts = np.arange(left_wing, right_wing + 1)
    # get probability for each integer within shifts
    probs = stats.norm.cdf(shifts + 0.5, mu, sigma) - stats.norm.cdf(shifts - 0.5, mu, sigma)
    # we need to add all the wings to the outermost bins
    probs[0] = stats.norm.cdf(left_wing + 0.5, mu, sigma)
    probs[-1] = 1 - stats.norm.cdf(right_wing - 0.5, mu, sigma)
    # check that sum of probabilities is 1
    assert (np.isclose(sum(probs), 1))
    probs *= amplitude_factor
    # put all the unused probability to bin with shift 0
    # so instead of rejecting the replenishment,
    # we will return shift 0, which will do nothing
    zero_index = np.where(shifts == 0)[0]
    probs[zero_index] += 1 - amplitude_factor
    # check that sum of probabilities is 1 again
    assert (np.isclose(sum(probs), 1))
    return lambda size: np.random.choice(shifts, p=probs, size=size)


def init_gen_0(nuc_num, init_den):
    gen_0_nuc_comp = np.zeros(nuc_num, dtype=int)
    for i in range(nuc_num):
        if random.random() < init_den:
            gen_0_nuc_comp[i] = 1
    return gen_0_nuc_comp


def replenish(old_nuc_comp, rr, mu, sigma, af):
    nuc_comp = np.copy(old_nuc_comp)
    for i in range(rr):
        random_shift = generate_probability_function(mu, sigma, af)
        cenp_a_positions = np.argwhere(nuc_comp == 1).flatten()
        distances_to_deposit = random_shift(cenp_a_positions.size)
        new_cenp_a_positions = cenp_a_positions + distances_to_deposit
        new_cenp_a_positions[new_cenp_a_positions >= len(nuc_comp)] -= len(nuc_comp)
        nuc_comp[new_cenp_a_positions] = 1
    return nuc_comp


def dilute(old_nuc_comp, distribution_bias):
    nuc_comp = np.copy(old_nuc_comp)
    cenp_a_positions = np.argwhere(nuc_comp == 1).flatten()
    for i in cenp_a_positions:
        if random.random() < distribution_bias:
            nuc_comp[i] = 0
    return nuc_comp


def initialize():
    global state
    state = init_gen_0(nucleosome_num, cenp_a_fraction)


def observe():
    global state
    density = state.sum() / len(state)
    positions = np.argwhere(state == 1).flatten()
    y = [1] * len(positions)
    plt.cla()
    plt.plot(positions, y, '.')
    plt.xlim(0, len(state))
    plt.title('density={}'.format(density))
    plt.axis('off')


def update():
    global state
    state = replenish(state, replenish_rounds, mu, sigma, amplitude_factor)
    state = dilute(state, dist_bias)


pycxsimulator.GUI(parameterSetters=[loading_efficiency, rounds_of_replenishment, loading_range, initial_density,
                                    number_of_nucleosomes, distribution_bias]).start(func=[initialize, observe, update])

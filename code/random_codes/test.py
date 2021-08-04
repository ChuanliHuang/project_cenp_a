import scipy.stats as stats
import numpy as np
import random
import matplotlib.pyplot as plt


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


def replenish(nuc_comp, _rs):
    for i in range(1):
        cenp_a_positions = np.argwhere(nuc_comp == 1).flatten()
        print(cenp_a_positions)
        distance_of_depositions = _rs(cenp_a_positions.size)
        print(distance_of_depositions)
        new_cenp_a_positions = cenp_a_positions + distance_of_depositions
        print(new_cenp_a_positions)
        new_cenp_a_positions[new_cenp_a_positions >= len(nuc_comp)] -= len(nuc_comp)
        print(new_cenp_a_positions)
        nuc_comp[new_cenp_a_positions] = 1
    return nuc_comp


def init_gen_0(nucleosome_num, fraction):
    gen_0_nuc_comp = np.zeros(nucleosome_num, dtype=int)
    for i in range(nucleosome_num):
        if random.random() < fraction:
            gen_0_nuc_comp[i] = 1
    return gen_0_nuc_comp

random_shift = generate_probability_function(0, 3, 1)
nuc_comp = init_gen_0(10, 0.5)
print(nuc_comp)
nuc_comp = replenish(nuc_comp, random_shift)
print(nuc_comp)


# plt.hist([int(random_shift(1)) for i in range(10000)], range=(-10,10), bins=21)
# plt.show()
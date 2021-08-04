import scipy.stats as stats
import numpy as np
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

rs = generate_probability_function(0, 3, 1)
res = rs(1000000)
plt.hist(res, bins=(res.max() - res.min()), density=True)
plt.show()

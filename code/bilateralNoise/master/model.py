import random
import numpy as np
from scipy import stats
import os


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


def replenish(old_state, rr, s, a, b):
    new_state = np.copy(old_state)
    f = generate_discrete_gaussian_distribution(s)
    for i in range(rr):
        for x in range(len(old_state)):
            if random.random() > b:
                p_incorp = 0
                for dx in range(3 * -s, 3 * s + 1):
                    p_incorp += old_state[((x + dx) % len(old_state))] * f[dx + 3 * s]
                p_incorp *= a
                if random.random() < p_incorp:
                    new_state[x] = 1
            else:
                if old_state[x] == 0:
                    new_state[x] = 1
                else:
                    new_state[x] = 0
        old_state = np.copy(new_state)
    return new_state


def dilute(old_state, bias):
    new_state = np.copy(old_state)
    one_positions = np.argwhere(old_state == 1).flatten()
    for i in one_positions:
        if random.random() < bias:
            new_state[i] = 0
    return new_state


def output_data(data, folder_name, lineage_num):
    # get current working directory
    cwd = os.getcwd()
    # create folder if not existing
    dir_name = cwd + '/' + folder_name
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    # create file name
    file_name = dir_name + '/lineage{}.csv'.format(lineage_num)
    # save as csv
    data = np.array(data)
    np.savetxt(file_name, data, delimiter=',')


if __name__ == '__main__':
    import time
    import yaml
    from bunch import bunchify

    start_time = time.perf_counter()
    with open("parameters.yaml") as f:
        config = yaml.load(f, Loader=yaml.Loader)
        config = bunchify(config)

    # simulate lineages
    for lineage in range(config.lineages_to_simulate):
        # initiate
        states = []
        state = init_gen_0(config.nucleosomeNumber, config.initialDensity)
        states.append(state)
        # simulate a lineage
        for generation in range(config.generations_to_simulate):
            # update
            state = replenish(state, config.replenish_rounds, config.sigma, config.alpha, config.beta)
            state = dilute(state, config.dilutionBias)
            # observe
            states.append(state)
            print('lineage {} generation {}'.format(lineage, generation))
        # output states
        output_data(states, 'states', lineage)

    end_time = time.perf_counter()
    rt = round(end_time - start_time, 2)
    print('Running time: {}s'.format(rt))
import scipy.stats as stats
import numpy as np
import random
import os


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


def replenish(old_nuc_comp, target_density, rr, mu, sigma, max_af):
    nuc_comp = np.copy(old_nuc_comp)
    af_list = []
    for i in range(rr):
        # detect current density
        density = nuc_comp.sum() / nuc_comp.shape[-1]
        # calculate af (linear relationship)
        af = max_af * (1 - density / target_density)
        if af < 0:
            af = 0
        af_list.append(af)
        random_shift = generate_probability_function(mu, sigma, af)
        cenp_a_positions = np.argwhere(nuc_comp == 1).flatten()
        distances_to_deposit = random_shift(cenp_a_positions.size)
        new_cenp_a_positions = cenp_a_positions + distances_to_deposit
        new_cenp_a_positions[new_cenp_a_positions >= len(nuc_comp)] -= len(nuc_comp)
        nuc_comp[new_cenp_a_positions] = 1
    return nuc_comp, af_list


def dilute(old_nuc_comp, distribution_bias):
    nuc_comp = np.copy(old_nuc_comp)
    cenp_a_positions = np.argwhere(nuc_comp == 1).flatten()
    for i in cenp_a_positions:
        if random.random() < distribution_bias:
            nuc_comp[i] = 0
    return nuc_comp


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
    from itertools import chain
    start_time = time.perf_counter()
    with open("parameters.yaml") as f:
        config = yaml.load(f, Loader=yaml.Loader)
        config = bunchify(config)

    # simulate lineages
    for lineage in range(config.simulation.lineages_to_simulate):
        # initiate
        states = []
        af_nums = []
        state = init_gen_0(config.simulation.lineage.nucleosome_num, config.simulation.lineage.cenp_a_fraction)
        states.append(state)
        # simulate a lineage
        for generation in range(config.simulation.lineage.generations_to_simulate):
            # update
            state, af_num = replenish(state,
                                      config.simulation.lineage.target_density,
                                      config.simulation.lineage.generation.replenish_rounds,
                                      config.simulation.lineage.generation.mu,
                                      config.simulation.lineage.generation.sigma,
                                      config.simulation.lineage.generation.amplitude_factor)
            state = dilute(state, config.simulation.lineage.generation.dist_bias)
            # observe

            # if (generation + 1) % 50 == 0:
            #     states.append(state)
            states.append(state)
            af_nums.append(af_num)
        # output states
        output_data(states, 'states', lineage)
        # output af_num
        af_nums = list(chain.from_iterable(af_nums))
        output_data(af_nums, 'af_num', lineage)
    end_time = time.perf_counter()
    rt = round(end_time - start_time, 2)
    print('Running time: {}s'.format(rt))

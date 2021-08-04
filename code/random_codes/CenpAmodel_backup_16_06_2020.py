import yaml
import numpy as np
import random
from bunch import bunchify
import scipy.stats as stats
import concurrent.futures
import copy
import time
import csv


class Simulation:
    def __init__(self, _config):
        self.config = _config

    def __repr__(self):
        return 'S={} RR={} AF={}'.format(self.config.simulation.lineage.generation.sigma,
                                         self.config.simulation.lineage.generation.replenishment_rounds,
                                         self.config.simulation.lineage.generation.amplitude_factor)

    def start(self, folder_name):
        random_shift = generate_probability_function(0, self.config.simulation.lineage.generation.sigma,
                                                     self.config.simulation.lineage.generation.amplitude_factor)
        # density_fname = '/home/chuang/project_cenp_a/data/' + folder_name + '/raw_data/density/' + str(
        #     self.config.simulation.lineage.generation.sigma) + '_' + str(
        #     self.config.simulation.lineage.generation.replenish_rounds) + '_' + str(
        #     self.config.simulation.lineage.generation.amplitude_factor) + '_' + str(
        #     self.config.simulation.lineage.cenp_a_fraction)
        density_fname = '/Users/kikawaryoku/PycharmProjects/project_cenp_a/replenish_efficiency/' + folder_name + '/raw_data/density/' + str(
            self.config.simulation.lineage.generation.sigma) + '_' + str(
            self.config.simulation.lineage.generation.replenish_rounds) + '_' + str(
            self.config.simulation.lineage.generation.amplitude_factor) + '_' + str(
            self.config.simulation.lineage.cenp_a_fraction)
        density_fname = density_fname.replace('.', '')
        density_fname = density_fname + '.csv'
        density_f = open(density_fname, 'a')
        density_wr = csv.writer(density_f, dialect='excel')
        # re_fname = '/home/chuang/project_cenp_a/data/' + folder_name + '/raw_data/replenish_efficiency/' + str(
        #     self.config.simulation.lineage.generation.sigma) + '_' + str(
        #     self.config.simulation.lineage.generation.replenish_rounds) + '_' + str(
        #     self.config.simulation.lineage.generation.amplitude_factor) + '_' + str(
        #     self.config.simulation.lineage.cenp_a_fraction)
        re_fname = '/Users/kikawaryoku/PycharmProjects/project_cenp_a/replenish_efficiency/' + folder_name + '/raw_data/replenish_efficiency/' + str(
            self.config.simulation.lineage.generation.sigma) + '_' + str(
            self.config.simulation.lineage.generation.replenish_rounds) + '_' + str(
            self.config.simulation.lineage.generation.amplitude_factor) + '_' + str(
            self.config.simulation.lineage.cenp_a_fraction)
        re_fname = re_fname.replace('.', '')
        re_fname = re_fname + '.csv'
        re_f = open(re_fname, 'a')
        re_wr = csv.writer(re_f, dialect='excel')
        for lineage_id in range(self.config.simulation.lineages_to_simulate):
            lineage = Lineage(self.config.simulation.lineage, lineage_id, random_shift)
            lineage.start(folder_name)
            re_wr.writerow(lineage.replenish_efficiencies)
            lineage_densities = []
            for _generation in lineage.generations:
                density = _generation.nuc_comp.sum() / self.config.simulation.lineage.nucleosome_num
                lineage_densities.append(density)
            density_wr.writerow(lineage_densities)
        re_f.close()
        density_f.close()


class Lineage:
    def __init__(self, _config, _id, _rs):
        self.config = _config
        self.id = _id
        self.rs = _rs
        self.generations = []
        self.replenish_efficiencies = []

    def __repr__(self):
        return 'Lineage {}'.format(self.id)

    def init_gen_0(self):
        gen_0_nuc_comp = np.zeros(self.config.nucleosome_num, dtype=int)
        for i in range(self.config.nucleosome_num):
            if random.random() < self.config.cenp_a_fraction:
                gen_0_nuc_comp[i] = 1
        gen_0 = Generation(self.config.generation, 0, gen_0_nuc_comp)
        return gen_0

    def start(self, folder_name):
        # file_name = '/home/chuang/project_cenp_a/data/' + folder_name + '/raw_data/nuc_comp/' + str(
        #     self.config.generation.sigma) + '_' + str(
        #     self.config.generation.replenish_rounds) + '_' + str(
        #     self.config.generation.amplitude_factor) + '_' + str(
        #     self.config.cenp_a_fraction) + '_' + str(self.id)
        file_name = '/Users/kikawaryoku/PycharmProjects/project_cenp_a/replenish_efficiency/' + folder_name + '/raw_data/nuc_comp/' + str(
            self.config.generation.sigma) + '_' + str(
            self.config.generation.replenish_rounds) + '_' + str(
            self.config.generation.amplitude_factor) + '_' + str(self.config.cenp_a_fraction) + '_' + str(self.id)
        file_name = file_name.replace('.', '')
        file_name = file_name + '.csv'
        lineage_f = open(file_name, 'a')
        lineage_wr = csv.writer(lineage_f, dialect='excel')
        gen = self.init_gen_0()
        self.generations.append(gen)
        lineage_wr.writerow(list(gen.nuc_comp))
        for gen_id in range(1, self.config.generations_to_simulate):
            nuc_comp, replenish_efficiency = gen.replicate(self.rs)
            self.replenish_efficiencies.append(replenish_efficiency)
            gen = Generation(self.config.generation, gen_id, nuc_comp)
            self.generations.append(gen)
            # if gen_id % 20 == 0:
            #     lineage_wr.writerow(list(gen.nuc_comp))
            lineage_wr.writerow(list(gen.nuc_comp))
        lineage_f.close()


class Generation:
    def __init__(self, _config, _id, _nuc_comp):
        self.config = _config
        self.id = _id
        self.nuc_comp = _nuc_comp

    def __repr__(self):
        return 'Generation {}\n'.format(self.id)

    def replenish(self, nuc_comp, _rs):
        for i in range(self.config.replenish_rounds):
            cenp_a_positions = np.argwhere(nuc_comp == 1).flatten()
            distances_to_deposit = _rs(cenp_a_positions.size)
            new_cenp_a_positions = cenp_a_positions + distances_to_deposit
            new_cenp_a_positions[new_cenp_a_positions >= len(nuc_comp)] -= len(nuc_comp)
            nuc_comp[new_cenp_a_positions] = 1
        return nuc_comp

    def dilute(self, nuc_comp):
        cenp_a_positions = np.argwhere(nuc_comp == 1).flatten()
        for i in cenp_a_positions:
            if random.random() < self.config.dist_bias:
                nuc_comp[i] = 0
        return nuc_comp

    def replicate(self, _rs):
        nuc_comp = np.copy(self.nuc_comp)
        old_sum = nuc_comp.sum()
        nuc_comp = self.replenish(nuc_comp, _rs)
        new_sum = nuc_comp.sum()
        replenish_efficiency = (new_sum - old_sum) / old_sum
        nuc_comp = self.dilute(nuc_comp)
        return nuc_comp, replenish_efficiency


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


def assign_range(_):
    if isinstance(_, tuple):
        _start, _stop, _step_size = _[0], _[1], _[2]
        _range = np.linspace(_start, _stop, _step_size)
    elif isinstance(_, list):
        _range = _
    else:
        _range = [_]
    return _range


def generate_parameter_combos(s, rr, af, cf):
    """kwargs can be tuple (start, stop, step), list or int"""
    s_range = assign_range(s)
    rr_range = assign_range(rr)
    af_range = assign_range(af)
    cf_range = assign_range(cf)
    parameter_combos = []
    for _s in s_range:
        for _rr in rr_range:
            for _af in af_range:
                for _cf in cf_range:
                    parameter_combo = (_s, _rr, _af, _cf)
                    parameter_combos.append(parameter_combo)
    return parameter_combos


def start_a_simulation(folder_name, parameter_combo):
    sigma, replenish_rounds, amplitude_factor, cenp_a_fraction = parameter_combo[0], parameter_combo[1], \
                                                                 parameter_combo[2], parameter_combo[3]
    _config = copy.copy(config)
    _config.simulation.lineage.generation.sigma = sigma
    _config.simulation.lineage.generation.replenish_rounds = replenish_rounds
    _config.simulation.lineage.generation.amplitude_factor = amplitude_factor
    _config.simulation.lineage.cenp_a_fraction = cenp_a_fraction
    simulation = Simulation(_config)
    simulation.start(folder_name)
    del simulation


def start(folder_name, s=3, rr=3, af=1, cf=0.02):
    parameter_combos = generate_parameter_combos(s, rr, af, cf)
    folder_names = [folder_name] * len(parameter_combos)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(start_a_simulation, folder_names, parameter_combos)


if __name__ == '__main__':
    start_time = time.perf_counter()
    with open("parameters.yaml") as f:
        config = yaml.load(f, Loader=yaml.Loader)
        config = bunchify(config)
    start('af=0.7_reversed', af=0.7, cf=(0.05, 1, 20))
    end_time = time.perf_counter()
    rt = round(end_time - start_time, 2)
    print('Running time: {}s'.format(rt))

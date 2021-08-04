import yaml
import numpy as np
import random
from bunch import bunchify
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
        file_name = '/home/chuang/project_cenp_a/data/' + folder_name + '/raw_data/density/' + str(
            self.config.simulation.lineage.generation.sigma) + '_' + str(
            self.config.simulation.lineage.generation.replenish_rounds) + '_' + str(
            self.config.simulation.lineage.generation.amplitude_factor)
        file_name = file_name.replace('.', '')
        file_name = file_name + '.csv'
        simulation_f = open(file_name, 'a')
        simulation_wr = csv.writer(simulation_f, dialect='excel')
        for lineage_id in range(self.config.simulation.lineages_to_simulate):
            lineage = Lineage(self.config.simulation.lineage, lineage_id)
            lineage.start(folder_name)
            lineage_densities = []
            for _generation in lineage.generations:
                density = _generation.nuc_comp.sum() / self.config.simulation.lineage.nucleosome_num
                lineage_densities.append(density)
            simulation_wr.writerow(lineage_densities)
        simulation_f.close()


class Lineage:
    def __init__(self, _config, _id):
        self.config = _config
        self.id = _id
        self.generations = []

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
        file_name = '/home/chuang/project_cenp_a/data/' + folder_name + '/raw_data/nuc_comp/' + str(
            self.config.generation.sigma) + '_' + str(
            self.config.generation.replenish_rounds) + '_' + str(
            self.config.generation.amplitude_factor) + '_' + str(self.id)
        file_name = file_name.replace('.', '')
        file_name = file_name + '.csv'
        lineage_f = open(file_name, 'a')
        lineage_wr = csv.writer(lineage_f, dialect='excel')
        gen = self.init_gen_0()
        self.generations.append(gen)
        lineage_wr.writerow(list(gen.nuc_comp))
        for gen_id in range(1, self.config.generations_to_simulate):
            nuc_comp = gen.replicate()
            gen = Generation(self.config.generation, gen_id, nuc_comp)
            self.generations.append(gen)
            if gen_id % 20 == 0:
                lineage_wr.writerow(list(gen.nuc_comp))
        lineage_f.close()


class Generation:
    def __init__(self, _config, _id, _nuc_comp):
        self.config = _config
        self.id = _id
        self.nuc_comp = _nuc_comp

    def __repr__(self):
        return 'Generation {}\n'.format(self.id)

    def replenish(self, nuc_comp):
        for i in range(self.config.replenish_rounds):
            cenp_a_positions = np.argwhere(nuc_comp == 1).flatten()
            for j in cenp_a_positions:
                # generate a random integer from a normal distribution as the distance of deposition
                rn = np.random.normal(self.config.mu, self.config.sigma)
                distance_of_deposition = round(rn)
                # deposit
                new_cenp_a_position = j + distance_of_deposition
                # avoid error raised because of boundary
                try:
                    if random.random() < self.config.amplitude_factor:
                        nuc_comp[new_cenp_a_position] = 1
                except IndexError:
                    pass
        return nuc_comp

    def dilute(self, nuc_comp):
        cenp_a_positions = np.argwhere(nuc_comp == 1).flatten()
        for i in cenp_a_positions:
            if random.random() < self.config.dist_bias:
                nuc_comp[i] = 0
        return nuc_comp

    def replicate(self):
        nuc_comp = np.copy(self.nuc_comp)
        nuc_comp = self.replenish(nuc_comp)
        nuc_comp = self.dilute(nuc_comp)
        return nuc_comp


def assign_range(_):
    if isinstance(_, tuple):
        start, stop, step_size = _[0], _[1], _[2]
        _range = np.linspace(start, stop, step_size)
    elif isinstance(_, list):
        _range = _
    else:
        _range = [_]
    return _range


def generate_parameter_combos(s=3, rr=3, af=1):
    '''kwargs can be tuple (start, stop, step), list or int'''
    S_range = assign_range(s)
    RR_range = assign_range(rr)
    AF_Range = assign_range(af)
    parameter_combos = []
    for _S in S_range:
        for _RR in RR_range:
            for _AF in AF_Range:
                parameter_combo = (_S, _RR, _AF)
                parameter_combos.append(parameter_combo)
    return parameter_combos


def start_a_simulation(folder_name, parameter_combo):
    sigma, replenish_rounds, amplitude_factor = parameter_combo[0], parameter_combo[1], parameter_combo[2]
    _config = copy.copy(config)
    _config.simulation.lineage.generation.sigma = sigma
    _config.simulation.lineage.generation.replenish_rounds = replenish_rounds
    _config.simulation.lineage.generation.amplitude_factor = amplitude_factor
    simulation = Simulation(_config)
    simulation.start(folder_name)
    del simulation


def Start(folder_name, s=3, rr=3, af=1):
    parameter_combos = generate_parameter_combos(s, rr, af)
    folder_names = [folder_name] * len(parameter_combos)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(start_a_simulation, folder_names, parameter_combos)


if __name__ == '__main__':
    start_time = time.perf_counter()
    with open("parameters.yaml") as f:
        config = yaml.load(f, Loader=yaml.Loader)
        config = bunchify(config)
    Start('exp_5', af=(0.5, 0.53, 11), rr=10)
    end_time = time.perf_counter()
    rt = round(end_time - start_time, 2)
    print('Running time: {}s'.format(rt))

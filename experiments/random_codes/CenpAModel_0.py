import yaml
import numpy as np
import random
from bunch import bunchify
import matplotlib.pyplot as plt
from itertools import chain
import concurrent.futures
import copy
import time
import pandas as pd


class Simulation:
    def __init__(self, _config):
        self.config = _config
        self.lineages = []

    def __repr__(self):
        return 'S={} RR={} AF={}'.format(self.config.simulation.lineage.generation.sigma,
                                         self.config.simulation.lineage.generation.replenishment_rounds,
                                         self.config.simulation.lineage.generation.amplitude_factor)

    def start(self):
        for lineage_id in range(self.config.simulation.lineages_to_simulate):
            lineage = Lineage(self.config.simulation.lineage, lineage_id)
            lineage.start()
            self.lineages.append(lineage)
            # print(lineage)
            # print('finished at {} {}:{}:{}'.format(time.localtime().tm_mday, time.localtime().tm_hour,
            #                                        time.localtime().tm_min, time.localtime().tm_sec))

    def output_density(self, folder_name):
        simulation_densities = []
        for _lineage in self.lineages:
            lineage_densities = []
            for _generation in _lineage.generations:
                density = _generation.nuc_comp.sum() / self.config.simulation.lineage.nucleosome_num
                lineage_densities.append(density)
            simulation_densities.append(lineage_densities)
        simulation_densities = pd.DataFrame(simulation_densities)
        file_name = '/home/chuang/project_cenp_a/data/' + folder_name + '/raw_data/density/' + str(
            self.config.simulation.lineage.generation.sigma) + '_' + str(
            self.config.simulation.lineage.generation.replenishment_rounds) + '_' + str(
            self.config.simulation.lineage.generation.amplitude_factor)
        file_name = file_name.replace('.', '')
        file_name = file_name + '.csv'
        simulation_densities.to_csv(file_name)
        del simulation_densities

    def show_density_trajectorys(self):
        for lineage in self.lineages:
            lineage.show_density_trajectory()
        plt.show()

    def show_kymographs(self):
        for lineage in self.lineages:
            lineage.show_kymograph()


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

    def start(self):
        gen = self.init_gen_0()
        self.generations.append(gen)
        for gen_id in range(1, self.config.generations_to_simulate):
            nuc_comp = gen.replicate()
            gen = Generation(self.config.generation, gen_id, nuc_comp)
            self.generations.append(gen)

    def output_nuc_comp(self, folder_name):
        compositions = []
        for _generation in self.generations:
            if _generation.id % 20 == 0:
                composition = list(_generation.nuc_comp)
                compositions.append(composition)
        compositions = pd.DataFrame(compositions)
        file_name = '/home/chuang/project_cenp_a/data/' + folder_name + '/raw_data/nuc_comp/' + str(
            self.config.generation.sigma) + '_' + str(
            self.config.generation.replenishment_rounds) + '_' + str(
            self.config.generation.amplitude_factor) + '_' + str(self.id)
        file_name = file_name.replace('.', '')
        file_name = file_name + '.csv'
        compositions.to_csv(file_name)
        del compositions

    def show_density_trajectory(self):
        x = range(self.config.generations_to_simulate)
        y = []
        for generation in self.generations:
            density = generation.nuc_comp.sum() / self.config.nucleosome_num
            y.append(density)
        plt.ylabel('density')
        plt.xlabel('generation')
        plt.plot(x, y)

    def show_kymograph(self):
        x = []
        y = []
        for generation in self.generations:
            x_coordinates = np.argwhere(generation.nuc_comp == 1).flatten().tolist()
            x.append(x_coordinates)
            y_coordinates = [generation.id] * generation.nuc_comp.sum()
            y.append(y_coordinates)
        x = list(chain.from_iterable(x))
        y = list(chain.from_iterable(y))
        plt.hist2d(x, y, bins=(self.config.nucleosome_num, self.config.generations_to_simulate), cmap='gray')
        plt.title('Lineage{}'.format(self.id))
        plt.ylabel('generation')
        plt.xlabel('nucleosome position')
        plt.show()


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
    '''kwargs can be tuple (start, stop, stop), list or int'''
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


def start_a_simulation(parameter_combo):
    sigma, replenishment_rounds, amplitude_factor = parameter_combo[0], parameter_combo[1], parameter_combo[2]
    _config = copy.copy(config)
    _config.simulation.lineage.generation.sigma = sigma
    _config.simulation.lineage.generation.replenishment_rounds = replenishment_rounds
    _config.simulation.lineage.generation.amplitude_factor = amplitude_factor
    simulation = Simulation(_config)
    simulation.start()
    simulation.show_density_trajectorys()
    simulation.show_kymographs()
    # simulation.output_density(folder_name)
    # for lineage in simulation.lineages:
        # lineage.output_nuc_comp(folder_name)
    del simulation


def Start(s=3, rr=3, af=1):
    parameter_combos = generate_parameter_combos(s, rr, af)
    # folder_names = [folder_name] * len(parameter_combos)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(start_a_simulation, parameter_combos)


if __name__ == '__main__':
    start_time = time.perf_counter()
    with open("parameters.yaml") as f:
        config = yaml.load(f, Loader=yaml.Loader)
        config = bunchify(config)
    Start(af=0.49)
    end_time = time.perf_counter()
    rt = round(end_time - start_time, 2)
    print('Running time: {}s'.format(rt))
    # s.show_kymographs()
    # s.show_density_trajectorys()

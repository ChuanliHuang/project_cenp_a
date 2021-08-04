import model
import time
import yaml
from bunch import bunchify
import os
from itertools import chain

start_time = time.perf_counter()
# change to the real folder
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)
with open("parameters.yaml") as f:
    config = yaml.load(f, Loader=yaml.Loader)
    config = bunchify(config)

# simulate lineages
for lineage in range(config.simulation.lineages_to_simulate):
    # initiate
    states = []
    af_nums = []
    state = model.init_gen_0(config.simulation.lineage.nucleosome_num, config.simulation.lineage.cenp_a_fraction)
    states.append(state)
    # simulate a lineage
    for generation in range(config.simulation.lineage.generations_to_simulate):
        # update
        state, af_num = model.replenish(state,
                                        config.simulation.lineage.target_density,
                                        config.simulation.lineage.generation.replenish_rounds,
                                        config.simulation.lineage.generation.mu,
                                        config.simulation.lineage.generation.sigma,
                                        config.simulation.lineage.generation.amplitude_factor)
        state = model.dilute(state, config.simulation.lineage.generation.dist_bias)
        # observe

        if (generation + 1) % 50 == 0:
            states.append(state)
        af_nums.append(af_num)
    # output states
    model.output_data(states, 'states', lineage)
    # output af_num
    af_nums = list(chain.from_iterable(af_nums))
    model.output_data(af_nums, 'af_num', lineage)
end_time = time.perf_counter()
rt = round(end_time - start_time, 2)
print('rr={} done'.format(config.simulation.lineage.generation.replenish_rounds))
print('Running time: {}s'.format(rt))

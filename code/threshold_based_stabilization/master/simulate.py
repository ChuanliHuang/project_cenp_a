import model
import time
import yaml
from bunch import bunchify
import os

start_time = time.perf_counter()
# change to the real folder
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)
with open("parameters.yaml") as f:
    config = yaml.load(f, Loader=yaml.Loader)
    config = bunchify(config)
# simulate lineages
for lineage in range(config.lineages_to_simulate):
    # initiate
    states = []
    state = model.init_gen_0(config.nucleosomeNumber, config.initialDensity)
    # save gen 0
    # states.append(state)
    # simulate a lineage
    for generation in range(config.generations_to_simulate):
        # update
        state = model.replenish(state, config.replenish_rounds, config.sigma, config.alpha, config.beta)
        state = model.dilute(state, config.dilutionBias)
        # observe
        if generation >= config.generations_to_simulate - 1000:  # only output the last 1000 generations
            states.append(state)
    # output states
    model.output_data(states, 'states', lineage)
end_time = time.perf_counter()
rt = round(end_time - start_time, 2)
print('alpha={} done'.format(config.alpha))
print('Running time: {}s'.format(rt))

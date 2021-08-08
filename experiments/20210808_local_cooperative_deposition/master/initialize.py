import numpy as np
import os
from pathlib import Path
from shutil import copy
import yaml


def figure_out_input_datatype(_):
    if isinstance(_, tuple):
        _start, _stop, _step_num = _[0], _[1], _[2]
        _range = np.linspace(_start, _stop, _step_num)
        _range = np.around(_range, 2)
    elif isinstance(_, list):
        _range = _
    else:
        _range = [_]
    return _range


def generate_parameter_combos(s, rr, af):
    """kwargs can be tuple (start, stop, step), list or num"""
    s_range = figure_out_input_datatype(s)
    rr_range = figure_out_input_datatype(rr)
    af_range = figure_out_input_datatype(af)
    parameter_combos = []
    for _s in s_range:
        for _rr in rr_range:
            for _af in af_range:
                parameter_combo = (_s, int(_rr), _af)
                parameter_combos.append(parameter_combo)
    return parameter_combos


def initialize_parameters_specific_folders(s=3, rr=3, af=1):
    """
    1. create parameters-specific folders in parent folder;
    2. copy model.py simulate.py parameters.yaml to each folder
    3. change parameters.yaml
    """
    parameter_combos = generate_parameter_combos(s, rr, af)
    # change to the real folder
    dir_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(dir_path)
    # get parent folder
    parent_path = Path(dir_path).parent
    # model.py pathname
    model_pathname = dir_path + '/' + 'model.py'
    # simulate.py pathname
    simulate_pathname = dir_path + '/' + 'simulate.py'
    # parameters.yaml pathname
    parameters_pathname = dir_path + '/' + 'parameters.yaml'
    # create parameters-specific folders
    for parameter_combo in parameter_combos:
        folder_name = 's{}rr{}af{}'.format(parameter_combo[0], parameter_combo[1], parameter_combo[2])
        folder_name = folder_name.replace('.', '')
        new_dir_path = str(parent_path) + '/results/' + folder_name
        if not os.path.exists(new_dir_path):
            os.makedirs(new_dir_path)
        # copy model.py simulate.py parameters.yaml to each folder
        copy(model_pathname, new_dir_path)
        copy(simulate_pathname, new_dir_path)
        copy(parameters_pathname, new_dir_path)
        # change parameters.yaml
        new_parameters_pathname = new_dir_path + '/' + 'parameters.yaml'
        with open(new_parameters_pathname) as f:
            list_doc = yaml.load(f, Loader=yaml.Loader)
            # change s
            list_doc['simulation']['lineage']['generation']['sigma'] = float(parameter_combo[0])
            # change rr
            list_doc['simulation']['lineage']['generation']['replenish_rounds'] = int(parameter_combo[1])
            # change af
            list_doc['simulation']['lineage']['generation']['amplitude_factor'] = float(parameter_combo[2])
        with open(new_parameters_pathname, "w") as f:
            yaml.dump(list_doc, f)
    print('parameters-specific folders initialized')


initialize_parameters_specific_folders(af=(0.85, 0.95, 11))

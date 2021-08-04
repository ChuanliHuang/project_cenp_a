import random
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from configurations import *
from tqdm import tqdm


def enter_a_cell(cell):
    """only enter cell with reasonable CENP-A density"""
    global data
    cenp_a_density = cell['nucleosome_composition'].sum() / nucleosome_num
    up_lim = (1 + density_max_variation) * cenp_a_fraction
    low_lim = (1 - density_max_variation) * cenp_a_fraction
    if low_lim < cenp_a_density < up_lim:
        data = np.append(data, [cell], axis=0)


def kill_cells(array):
    up_lim = (1 + density_max_variation) * cenp_a_fraction
    low_lim = (1 - density_max_variation) * cenp_a_fraction
    density_ls = []
    for i in array:
        cenp_a_density = i['nucleosome_composition'].sum() / nucleosome_num
        density_ls.append(cenp_a_density)
    density_array = np.array(density_ls)
    selection = ((density_array >= low_lim) & (density_array <= up_lim))
    viable_cell_positions = np.argwhere(selection).flatten()
    array = array[viable_cell_positions]
    return array


def init_cell_0():
    """initiate a certain number of cell 0"""
    global data
    for i in range(initiate_cell_number):
        cell_0_nuc_comp = np.zeros(nucleosome_num, dtype=int)
        for j in range(nucleosome_num):
            if random.random() < cenp_a_fraction:
                cell_0_nuc_comp[j] = 1
        cell_0 = np.array((i, 0, 0, cell_0_nuc_comp), dtype=dt)
        data = np.append(data, [cell_0], axis=0)
    data = kill_cells(data)


def replenish_cenp_a(mother):
    """replenishment that happens in a positional-based method"""
    cenp_a_positions = np.argwhere(mother['nucleosome_composition'] == 1).flatten()

    for i in cenp_a_positions:
        # generate a random integer from a normal distribution as the distance of deposition
        rn = np.random.normal(mu, sigma)
        distance_of_deposition = round(rn)
        # deposit
        new_cenp_a_position = i + distance_of_deposition
        # if new_cenp_a_position in range(nucleosome_num):
        try:
            if random.random() < amplitude_factor:
                mother['nucleosome_composition'][new_cenp_a_position] = 1
        except IndexError:
            pass
    return mother


def dilute_cenp_a(mother):
    """CENP-A from mother strand distributed to daughter strands"""
    dgt_1_nuc_comp = np.zeros(nucleosome_num, dtype=int)
    dgt_2_nuc_comp = np.zeros(nucleosome_num, dtype=int)
    for i in range(nucleosome_num):
        if mother['nucleosome_composition'][i] == 1:
            if random.random() < dist_bias:
                dgt_1_nuc_comp[i] = 1
            else:
                dgt_2_nuc_comp[i] = 1
    return dgt_1_nuc_comp, dgt_2_nuc_comp


def replicate_a_cell(cell):
    """replicate a selected cell"""
    global data
    mother = cell.copy()
    dgt_generation = mother['generation'] + 1
    dgt_parent = cell['cell_id']
    # replenishment
    for i in range(replenish_rounds):
        mother = replenish_cenp_a(mother)
    # dilution
    dgt_1_nuc_comp, dgt_2_nuc_comp = dilute_cenp_a(mother)
    dgt1_id = initiate_cell_number + 2 * dgt_parent
    dgt2_id = initiate_cell_number + 2 * dgt_parent + 1
    dgt1 = np.array((dgt1_id, dgt_generation, dgt_parent, dgt_1_nuc_comp), dtype=dt)
    dgt2 = np.array((dgt2_id, dgt_generation, dgt_parent, dgt_2_nuc_comp), dtype=dt)
    return dgt1, dgt2


def replicate_a_generation(generation_num):
    global data
    gen_data = np.array([], dtype=dt)
    for i in tqdm(data[data['generation'] == generation_num]):
        dgt1, dgt2 = replicate_a_cell(i)
        gen_data = np.append(gen_data, [dgt1], axis=0)
        gen_data = np.append(gen_data, [dgt2], axis=0)
    # plot_gen_distribution(gen_data)
    gen_data = kill_cells(gen_data)
    data = np.append(data, gen_data)


def plot_gen_distribution(gen_data):
    density_ls = []
    for j in gen_data:
        cenp_a_density = j['nucleosome_composition'].sum() / nucleosome_num
        density_ls.append(cenp_a_density)
    ax = sns.distplot(density_ls, kde=False)
    plt.title('AF={} RR={} S={}'.format(amplitude_factor, replenish_rounds, sigma))
    plt.show()
    plt.clf()


def plot_gen_average():
    global data
    for i in range(generations_to_simulate + 1):
        gen_data = data[data['generation'] == i]
        density_mean = gen_data['nucleosome_composition'].sum() / (nucleosome_num * len(gen_data))
        plt.plot(i, density_mean, color='b', marker='.')
    plt.xlabel('generation')
    plt.ylabel('CENP_A density_mean')
    plt.title('AF={} RR={} S={}'.format(amplitude_factor, replenish_rounds, sigma))
    plt.show()


def start():
    """start simulation"""
    global data
    init_cell_0()
    for j in range(generations_to_simulate):
        replicate_a_generation(j)


def show_APEX_plot():
    plt.imshow(data['nucleosome_composition'], cmap='gray')
    plt.title('AF={} RR={} MV={}'.format(amplitude_factor, replenish_rounds,
                                         density_max_variation))
    plt.ylabel('cell id')
    plt.xlabel('nucleosome position')
    white_patch = mpatches.Patch(color='white', label='CENP-A')
    black_patch = mpatches.Patch(color='black', label='H3')
    plt.legend(handles=[white_patch, black_patch])
    file_name = r'C:\Users\s1943350\PycharmProjects\CENP_A_modelling\figures\exp_2_replenishment\5\APEX_like.png'
    plt.savefig(file_name)
    plt.clf()


def show_cenp_a_density():
    sns.set(style="whitegrid")
    density_ls = []
    generation_ls = data['generation'].tolist()
    for i in range(len(data)):
        density = data['nucleosome_composition'][i].sum() / nucleosome_num
        density_ls.append(density)
    d = {'generation': generation_ls, 'CENP-A density': density_ls}
    df = pd.DataFrame(d)
    ax = sns.swarmplot(x="generation", y="CENP-A density", data=df)
    plt.title('S={} AF={} RR={}'.format(sigma, amplitude_factor, replenish_rounds))
    plt.show()
    # file_name = r'C:\Users\s1943350\PycharmProjects\CENP_A_modelling\figures\exp_2_replenishment\5\cenp_a_density.png'
    # plt.savefig(file_name)
    # plt.clf()


def save_data():
    df = pd.DataFrame({'generation': data['generation'], 'parent': data['parent'],
                       'nucleosome_composition': data['nucleosome_composition'].tolist()})
    df.index.name = 'cell_id'
    file_name = r'C:\Users\s1943350\PycharmProjects\CENP_A_modelling\figures\exp_2_replenishment\5\raw_data.xlsx'
    df.to_excel(file_name)


def save_configurations():
    with open(r'C:\Users\s1943350\PycharmProjects\CENP_A_modelling\configurations.py') as f:
        with open(
                r'C:\Users\s1943350\PycharmProjects\CENP_A_modelling\figures\exp_2_replenishment\5\configurations.txt',
                'w') as f1:
            for line in f:
                f1.write(line)


if __name__ == '__main__':
    start()
    plot_gen_average()
    # show_cenp_a_density()
    # show_APEX_plot()
    # save_data()
    # save_configurations()

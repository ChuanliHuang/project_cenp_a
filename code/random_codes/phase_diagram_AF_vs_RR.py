import random
import matplotlib.pyplot as plt
from configurations_AF_vs_RR import *
from tqdm import tqdm
import seaborn as sns


def init_cell_0():
    """initiate a certain number of cell 0"""
    cell_0_nuc_comp = np.zeros(nucleosome_num, dtype=int)
    for i in range(nucleosome_num):
        if random.random() < cenp_a_fraction:
            cell_0_nuc_comp[i] = 1
    cell_0 = np.array((0, 0, cell_0_nuc_comp), dtype=dt_1_dgt)
    return cell_0


def replenish_cenp_a(mother):
    """use previous results as template"""
    for j in range(replenish_rounds):
        cenp_a_positions = np.argwhere(mother['nucleosome_composition'] == 1).flatten()
        for i in cenp_a_positions:
            # generate a random integer from a normal distribution as the distance of deposition
            rn = np.random.normal(mu, sigma)
            distance_of_deposition = round(rn)
            # deposit
            new_cenp_a_position = i + distance_of_deposition
            try:
                if random.random() < amplitude_factor:
                    mother['nucleosome_composition'][new_cenp_a_position] = 1
            except IndexError:
                pass
    return mother



def dilute_cenp_a(mother):
    """CENP-A from mother strand distributed to daughter strands"""
    dgt_nuc_comp = np.zeros(nucleosome_num, dtype=int)
    for i in range(nucleosome_num):
        if mother['nucleosome_composition'][i] == 1:
            if random.random() < dist_bias:
                dgt_nuc_comp[i] = 1
    return dgt_nuc_comp


def replicate_a_cell(mother):
    """replicate a selected cell"""
    dgt_generation = mother['generation'] + 1
    lineage = mother['lineage']
    # replenishment
    mother = replenish_cenp_a(mother)
    # dilution
    dgt_nuc_comp = dilute_cenp_a(mother)
    dgt = np.array((dgt_generation, lineage, dgt_nuc_comp), dtype=dt_1_dgt)
    return dgt


def start():
    global data
    data = np.array([], dtype=dt_1_dgt)  # an array to store structured arrays
    cell = init_cell_0()
    data = np.append(data, [cell], axis=0)
    for i in range(generations_to_simulate):
        cell = replicate_a_cell(cell)
        data = np.append(data, [cell], axis=0)




if __name__ == '__main__':
    # arr = np.empty([af_step, RR_step])
    # for i in tqdm(range(af_step)):
    #     amplitude_factor = np.linspace(af_start, af_stop, af_step)[i]
    #     for j in range(RR_step):
    #         replenish_rounds = int(np.linspace(RR_start, RR_stop, RR_step)[j])
    #         start()
    #         first_hundred_sum = data[0:100]['nucleosome_composition'].sum()
    #         last_hundred_sum = data[generations_to_simulate - 100:generations_to_simulate]['nucleosome_composition'].sum()
    #         if last_hundred_sum >= first_hundred_sum:
    #             arr[i][j] = 1
    #         else:
    #             arr[i][j] = 0
    # plt.imshow(arr, cmap='gray')
    # # plt.colorbar()
    # x_ls = np.linspace(RR_start, RR_stop, 11).round(2).tolist()
    # y_ls = np.linspace(af_start, af_stop, 11).round(2).tolist()
    # x_ticks = map(str, x_ls)
    # y_ticks = map(str, y_ls)
    # plt.xlabel('replenishment rounds')
    # plt.xticks(np.linspace(0, RR_step - 1, 11), x_ticks)
    # plt.yticks(np.linspace(0, af_step - 1, 11), y_ticks)
    # plt.ylabel('loading efficiency')
    # plt.show()
    tran_AFs = []
    first_stabilized_densitys = []
    for j in tqdm(range(RR_step)):
        replenish_rounds = int(np.linspace(RR_start, RR_stop, RR_step)[j])
        for i in range(af_step):
            amplitude_factor = np.linspace(af_start, af_stop, af_step)[i]
            start()
            first_hundred_sum = data[0:100]['nucleosome_composition'].sum()
            last_hundred_sum = data[generations_to_simulate - 100:generations_to_simulate]['nucleosome_composition'].sum()
            if last_hundred_sum >= first_hundred_sum:
                tran_AFs.append(amplitude_factor)
                density = data[generations_to_simulate - 200:generations_to_simulate][
                                   'nucleosome_composition'].sum() / (200 * nucleosome_num)
                first_stabilized_densitys.append(density)
                break
    # create figure and axis objects with subplots()
    fig, ax = plt.subplots()
    # make a plot
    ax.plot(np.linspace(RR_start, RR_stop, RR_step), tran_AFs, color="red", marker="o")
    # set x-axis label
    ax.set_xlabel('replenishment rounds')
    # set y-axis label
    ax.set_ylabel("transition AF", color="red")
    # # twin object for two different y-axis on the sample plot
    # ax2 = ax.twinx()
    # # make a plot with different y-axis using second axis object
    # ax2.plot(np.linspace(RR_start, RR_stop, RR_step), first_stabilized_densitys, color="blue", marker="o")
    # ax2.set_ylabel("transition CENP-A density", color="blue")
    plt.title('Sigma = {}'.format(sigma))
    print(tran_AFs)
    plt.show()



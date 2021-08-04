import random
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from configurations_1_dgt import *
from tqdm import tqdm
from itertools import chain


def init_cell_0(j):
    """initiate a certain number of cell 0"""
    cell_0_nuc_comp = np.zeros(nucleosome_num, dtype=int)
    for i in range(nucleosome_num):
        if random.random() < cenp_a_fraction:
            cell_0_nuc_comp[i] = 1
    cell_0 = np.array((0, j, cell_0_nuc_comp), dtype=dt_1_dgt)
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


def replenish_cenp_a_new(mother):
    """must fulfill flanking"""
    cenp_a_positions = np.argwhere(mother['nucleosome_composition'] == 1).flatten()
    for i in cenp_a_positions:
        # local deposition
        if random.random() < amplitude_factor:
            try:
                mother['nucleosome_composition'][i + 1] = 1
            except:
                pass
            else:
                mother['nucleosome_composition'][i - 1] = 1
            finally:
                pass
        # long range deposition
        if random.random() < amplitude_factor_long:
            try:
                mother['nucleosome_composition'][random.randrange(0, nucleosome_num)] = 1
            except:
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
    for j in tqdm(range(initiate_cell_number)):
        cell = init_cell_0(j)
        data = np.append(data, [cell], axis=0)
        for i in range(generations_to_simulate):
            cell = replicate_a_cell(cell)
            data = np.append(data, [cell], axis=0)
    return data


def show_density_trajectory():
    for i in range(initiate_cell_number):
        lineage = data[data['lineage'] == i]
        x = []
        y = []
        for j in lineage:
            generation = int(j['generation'])
            density = j['nucleosome_composition'].sum() / nucleosome_num
            x.append(generation)
            y.append(density)
        plt.plot(x, y)
    plt.xlabel('generation')
    plt.ylabel('CENP-A density')
    plt.show()
    # file_name = '/Users/kikawaryoku/PycharmProjects/project_cenp_a/exp_3/density_evolution/' + str(trial) + '.png'
    # plt.savefig(file_name)
    # plt.clf()


def save_data():
    df = pd.DataFrame({'generation': data['generation'], 'lineage': data['lineage'],
                       'nucleosome_composition': data['nucleosome_composition'].tolist()})
    df.index.name = 'cell_id'
    file_name = '/Users/kikawaryoku/PycharmProjects/project_cenp_a/exp_3/raw_data/' + str(trial) + '.csv'
    df.to_csv(file_name)


def save_configurations():
    file_name = '/Users/kikawaryoku/PycharmProjects/project_cenp_a/exp_3/configurations/' + str(trial) + '.txt'
    with open('/random_codes/configurations_1_dgt.py') as f:
        with open(file_name, 'w') as f1:
            for line in f:
                f1.write(line)


def show_APEX_plot():
    plt.imshow(data['nucleosome_composition'], cmap='gray')
    plt.axis('off')
    plt.xlabel('nucleosome position')
    # white_patch = mpatches.Patch(color='white', label='CENP-A')
    # black_patch = mpatches.Patch(color='black', label='H3')
    # plt.legend(handles=[white_patch, black_patch])
    file_name = '/Users/kikawaryoku/PycharmProjects/project_cenp_a/exp_3/APEX_like_plots/' + str(trial) + '.png'
    plt.savefig(file_name)
    plt.show()
    plt.clf()


def show_example_APEX():
    cell_0 = data[0]['nucleosome_composition']
    cell_1 = data[int(generations_to_simulate / 4)]['nucleosome_composition']
    cell_2 = data[int(generations_to_simulate / 2)]['nucleosome_composition']
    cell_3 = data[int(generations_to_simulate * 3 / 4)]['nucleosome_composition']
    cell_4 = data[generations_to_simulate]['nucleosome_composition']
    ls = [cell_0, cell_1, cell_2, cell_3, cell_4]
    arr = np.array(ls)
    arr = np.repeat(arr, int(generations_to_simulate / 5), axis=0)
    plt.imshow(arr, cmap='gray')
    plt.axis('off')
    plt.xlabel('nucleosome position')
    plt.show()


def show_hist2d():
    x = []
    for i in data:
        arr = np.argwhere(i['nucleosome_composition'] == 1).flatten()
        arr = arr.tolist()
        x.append(arr)
    x = list(chain.from_iterable(x))
    y = []
    for j in range(generations_to_simulate + 1):
        n = data[j]['nucleosome_composition'].sum()
        ls = [j] * n
        y.append(ls)
    y = list(chain.from_iterable(y))
    plt.hist2d(x, y, bins=(nucleosome_num, generations_to_simulate), cmap='gray')
    plt.title('AF = {} S = {}'.format(amplitude_factor, sigma))
    plt.ylabel('generation')
    plt.xlabel('nucleosome position')
    plt.show()


if __name__ == '__main__':
    # start()
    # last_200_means = []
    # for i in range(initiate_cell_number):
    #     lineage = data[data['lineage'] == i]
    #     last_200_sum = lineage[generations_to_simulate - 200:generations_to_simulate]['nucleosome_composition'].sum()
    #     last_200_mean = last_200_sum / (200 * nucleosome_num)
    #     last_200_means.append(last_200_mean)
    #     first_hundred_sum = lineage[0:100]['nucleosome_composition'].sum()
    #     last_hundred_sum = lineage[generations_to_simulate - 100:generations_to_simulate]['nucleosome_composition'].sum()
    #     x = []
    #     y = []
    #     for j in lineage:
    #         generation = int(j['generation'])
    #         density = j['nucleosome_composition'].sum() / nucleosome_num
    #         x.append(generation)
    #         y.append(density)

    # x = range(generations_to_simulate + 1)
    # y = []
    # for j in range(generations_to_simulate + 1):
    #     generation = data[data['generation'] == j]
    #     density_mean = generation['nucleosome_composition'].sum() / (nucleosome_num * initiate_cell_number)
    #     y.append(density_mean)
    # plt.plot(x, y)
    # # final_mean = round(sum(last_200_means) / initiate_cell_number, 2)
    # plt.title('AF={} RR={} S={}'.format(amplitude_factor, replenish_rounds, sigma))
    # plt.xlabel('generation')
    # plt.ylabel('CENP-A density')
    # plt.show()

    # show_hist2d()
    # save_data()
    # save_configurations()
    # show_density_trajectory()
    # show_APEX_plot()
    # sigmas_of_interest = [1.4, 2.2, 3]
    # AF_resolution = 11
    # AF_range_0 = np.linspace(0.81, 0.82, AF_resolution)
    # AF_range_1 = np.linspace(0.6, 0.61, AF_resolution)
    # AF_range_2 = np.linspace(0.51, 0.52, AF_resolution)
    # AF_ranges = [AF_range_0, AF_range_1, AF_range_2]
    # final_df = pd.DataFrame(columns=['RR', 'S', 'AF', 'lineage', 'stabilized_1000', 'stabilized_2000', 'mean_1000', 'mean_2000'])
    # for i in range(len(sigmas_of_interest)):
    #     sigma = sigmas_of_interest[i]
    #     AF_range = AF_ranges[i]
    #     for amplitude_factor in tqdm(AF_range):
    #         data = np.array([], dtype=dt_1_dgt)
    #         data = start()
    #         for j in range(initiate_cell_number):
    #             lineage = data[data['lineage'] == j]
    #             first_hundred_sum = lineage[0:100]['nucleosome_composition'].sum()
    #             last_hundred_sum_1000 = lineage[900:1000]['nucleosome_composition'].sum()
    #             last_hundred_sum_2000 = lineage[1900:2000]['nucleosome_composition'].sum()
    #             decision_1000 = last_hundred_sum_1000 >= first_hundred_sum
    #             decision_2000 = last_hundred_sum_2000 >= first_hundred_sum
    #             mean_1000 = lineage[800:1000]['nucleosome_composition'].sum() / (200 * nucleosome_num)
    #             mean_2000 = lineage[1800:2000]['nucleosome_composition'].sum() / (200 * nucleosome_num)
    #             lineage_dict = {
    #                 'RR': replenish_rounds,
    #                 'S': sigma,
    #                 'AF': amplitude_factor,
    #                 'lineage': j,
    #                 'stabilized_1000':decision_1000,
    #                 'stabilized_2000':decision_2000,
    #                 'mean_1000':mean_1000,
    #                 'mean_2000':mean_2000
    #             }
    #             lineage_df = pd.DataFrame(lineage_dict, index=[0])
    #             final_df = final_df.append(lineage_dict, ignore_index=True)
    # file_name = '/Users/kikawaryoku/PycharmProjects/project_cenp_a/exp_4/raw_data/big_table.csv'
    # final_df.to_csv(file_name)

    # sigmas_of_interest = [1.4, 2.2, 3]
    # AF_resolution = 3
    # AF_range_0 = np.linspace(0.81, 0.82, AF_resolution)
    # AF_range_1 = np.linspace(0.6, 0.61, AF_resolution)
    # AF_range_2 = np.linspace(0.51, 0.52, AF_resolution)
    # AF_ranges = [AF_range_0, AF_range_1, AF_range_2]
    # for i in range(len(sigmas_of_interest)):
    #     sigma = sigmas_of_interest[i]
    #     AF_range = AF_ranges[i]
    #     for amplitude_factor in tqdm(AF_range):
    #         data = np.array([], dtype=dt_1_dgt)
    #         data = start()
    #         for j in range(initiate_cell_number):
    #             lineage = data[data['lineage'] == j]
    #             for k in lineage['nucleosome_composition']:
    #                 density = k.sum() / nucleosome_num

    # sigmas_of_interest = [1.4, 2.2, 3]
    # AF_resolution = 11
    # AF_range_0 = np.linspace(0.81, 0.82, AF_resolution)
    # AF_range_1 = np.linspace(0.6, 0.61, AF_resolution)
    # AF_range_2 = np.linspace(0.51, 0.52, AF_resolution)
    # AF_ranges = [AF_range_0, AF_range_1, AF_range_2]
    # cm = plt.get_cmap('gist_rainbow')
    # for k in range(len(sigmas_of_interest)):
    #     sigma = sigmas_of_interest[k]
    #     AF_range = AF_ranges[k]
    #     for i in tqdm(range(AF_resolution)):
    #         amplitude_factor = AF_range[i]
    #         data = np.array([], dtype=dt_1_dgt)
    #         data = start()
    #         x = range(generations_to_simulate + 1)
    #         y = []
    #         for j in range(generations_to_simulate + 1):
    #             generation = data[data['generation'] == j]
    #             density_mean = generation['nucleosome_composition'].sum() / (nucleosome_num * initiate_cell_number)
    #             y.append(density_mean)
    #         z = np.polyfit(x, y, 1)
    #         p = np.poly1d(z)
    #         # plt.scatter(x, y, color=cm(1.*i/AF_resolution), marker='.')
    #         plt.plot(x, p(x), color=cm(1.*i/AF_resolution), label='AF={} a={}'.format(round(amplitude_factor, 3), z[0]))
    #     plt.title('S={}'.format(sigma))
    #     plt.ylabel('CENP-A density')
    #     plt.xlabel('generation')
    #     plt.legend(prop={"size":7})
    #     file_name = '/home/chuang/project_cenp_a/data/{}.png'.format(k)
    #     plt.savefig(file_name, dpi=200)
    #     plt.clf()

    # sigmas_of_interest = [1.4, 2.2, 3]
    # AF_resolution = 2
    # AF_range_0 = np.linspace(0.81, 0.82, AF_resolution)
    # AF_range_1 = np.linspace(0.6, 0.61, AF_resolution)
    # AF_range_2 = np.linspace(0.51, 0.52, AF_resolution)
    # AF_ranges = [AF_range_0, AF_range_1, AF_range_2]
    # cm = plt.get_cmap('gist_rainbow')
    # for k in range(len(sigmas_of_interest)):
    #     sigma = sigmas_of_interest[k]
    #     AF_range = AF_ranges[k]
    #     for i in tqdm(range(AF_resolution)):
    #         amplitude_factor = AF_range[i]
    #         data = np.array([], dtype=dt_1_dgt)
    #         data = start()
    #         x = range(generations_to_simulate + 1)
    #         y = []
    #         for j in range(generations_to_simulate + 1):
    #             generation = data[data['generation'] == j]
    #             density_mean = generation['nucleosome_composition'].sum() / (nucleosome_num * initiate_cell_number)
    #             y.append(density_mean)
    #         plt.plot(x, y, color=cm(1.*i/AF_resolution), label='AF={}'.format(round(amplitude_factor, 3)))
    #     plt.title('S={}'.format(sigma))
    #     plt.ylabel('CENP-A density')
    #     plt.xlabel('generation')
    #     plt.legend()
    #     plt.show()
    #     file_name = '/home/chuang/project_cenp_a/data/{}.png'.format(k)
    #     plt.savefig(file_name, dpi=200)
    #     plt.clf()

    # sigmas_of_interest = [1.4, 2.2, 3]
    # AF_resolution = 11
    # AF_range_0 = np.linspace(0.82, 1, AF_resolution)
    # AF_range_1 = np.linspace(0.61, 1, AF_resolution)
    # AF_range_2 = np.linspace(0.52, 1, AF_resolution)
    # AF_ranges = [AF_range_0, AF_range_1, AF_range_2]
    # cm = plt.get_cmap('gist_rainbow')
    # for k in range(len(sigmas_of_interest)):
    #     sigma = sigmas_of_interest[k]
    #     AF_range = AF_ranges[k]
    #     for i in tqdm(range(AF_resolution)):
    #         amplitude_factor = AF_range[i]
    #         data = np.array([], dtype=dt_1_dgt)
    #         data = start()
    #         x = range(generations_to_simulate + 1)
    #         y = []
    #         for j in range(generations_to_simulate + 1):
    #             generation = data[data['generation'] == j]
    #             density_mean = generation['nucleosome_composition'].sum() / (nucleosome_num * initiate_cell_number)
    #             y.append(density_mean)
    #         plt.plot(x, y, color=cm(1.*i/AF_resolution), label='AF={}'.format(round(amplitude_factor, 3)))
    #     plt.title('S={}'.format(sigma))
    #     plt.ylabel('CENP-A density')
    #     plt.xlabel('generation')
    #     plt.legend()
    #     file_name = '/Users/kikawaryoku/PycharmProjects/project_cenp_a/exp_4/raw_data/when_stabilized_{}.png'.format(k)
    #     plt.savefig(file_name, dpi=200)
    #     plt.clf()

    # sigmas_of_interest = [1.4, 2.2, 3]
    # AF_resolution = 20
    # AF_range_0 = np.linspace(0.82, 1, AF_resolution)
    # AF_range_1 = np.linspace(0.61, 1, AF_resolution)
    # AF_range_2 = np.linspace(0.52, 1, AF_resolution)
    # AF_ranges = [AF_range_0, AF_range_1, AF_range_2]
    # for k in range(len(sigmas_of_interest)):
    #     sigma = sigmas_of_interest[k]
    #     AF_range = AF_ranges[k]
    #     last_250_means = []
    #     for i in tqdm(range(AF_resolution)):
    #         amplitude_factor = AF_range[i]
    #         data = np.array([], dtype=dt_1_dgt)
    #         data = start()
    #         generation_means = []
    #         for j in range(generations_to_simulate + 1):
    #             generation = data[data['generation'] == j]
    #             generation_mean = generation['nucleosome_composition'].sum() / (nucleosome_num * initiate_cell_number)
    #             generation_means.append(generation_mean)
    #         last_250_mean = sum(generation_means[generations_to_simulate-250:generations_to_simulate]) / 250
    #         last_250_means.append(last_250_mean)
    #     plt.scatter(AF_range, last_250_means)
    #     plt.title('S={}'.format(sigma))
    #     plt.ylabel('CENP-A density')
    #     plt.xlabel('loading efficiency')
    #     file_name = '/home/chuang/project_cenp_a/data/AF_test_{}.png'.format(k)
    #     plt.savefig(file_name, dpi=200)
    #     plt.clf()
    # data = start()
    # show_density_trajectory()
    # generation_means = []
    # for j in range(generations_to_simulate + 1):
    #     generation = data[data['generation'] == j]
    #     generation_mean = generation['nucleosome_composition'].sum() / (nucleosome_num * initiate_cell_number)
    #     generation_means.append(generation_mean)
    # last_250_mean = sum(generation_means[generations_to_simulate-250:generations_to_simulate]) / 250
    # print(last_250_mean)
    sigma = 1.4
    AF_resolution = 51
    AF_range = np.linspace(0.5, 1, AF_resolution)
    last_2000_means =[]
    for i in tqdm(range(AF_resolution)):
        amplitude_factor = AF_range[i]
        data = np.array([], dtype=dt_1_dgt)
        data = start()
        generation_means = []
        for j in range(generations_to_simulate + 1):
            generation = data[data['generation'] == j]
            generation_mean = generation['nucleosome_composition'].sum() / (nucleosome_num * initiate_cell_number)
            generation_means.append(generation_mean)
        last_2000_mean = sum(generation_means[generations_to_simulate-2000:generations_to_simulate]) / 2000
        last_2000_means.append(last_2000_mean)
    plt.scatter(AF_range, last_2000_means)
    plt.title('S={} RR={}'.format(sigma, replenish_rounds))
    plt.ylabel('CENP-A density')
    plt.xlabel('loading efficiency')
    file_name = '/home/chuang/project_cenp_a/data/AF_vs_density_0.png'
    plt.savefig(file_name, dpi=200)
    plt.clf()
import pandas as pd
from scipy.stats import pearsonr
from utils.codon_aa_dict import CODON_SIZE, AA_TO_CODON_DICT, CODON_TO_AA_DICT


def delta_min_max(minmax1, minmax2):
    """
    :param minmax1: minmax profile of protein 1 (a list of values)
    :param minmax2: minmax profile of protein 2 (a list of values)
    :return: the sum of the  absolute difference between each value of both lists divided by the number of values
    """
    if len(minmax1) != len(minmax2):
        raise ValueError('Both MinMax profiles must have the same length')

    absolute_differences = []
    for i in range(len(minmax1)):
        absolute_differences.append(abs(minmax1[i]-minmax2[i]))

    return sum(absolute_differences)/len(absolute_differences)


def correlation_min_max(minmax1, minmax2):
    """
    :param minmax1: minmax profile of protein 1 (a list of values)
    :param minmax2: minmax profile of protein 2 (a list of values)
    :return: the pearson correlation between minmax1 and minmax2 values
    """
    if len(minmax1) != len(minmax2):
        raise ValueError('Both MinMax profiles must have the same length')

    list1 = pd.Series(minmax1)
    list2 = pd.Series(minmax2)

    corr, _ = pearsonr(list1, list2)
    return corr


def average_min_max_codon_freq_per_aa(cu_table):
    average_codon_freq_per_aa_dict = dict()
    for aa, codons in AA_TO_CODON_DICT.items():
        codon_values = [codon_val for codon, codon_val in cu_table.get_codon_usage_table().items() if
                        codon in codons]
        aa_avg_min_max = dict()
        aa_avg_min_max['average'] = sum(codon_values) / len(codon_values)
        aa_avg_min_max['min'] = min(codon_values)
        aa_avg_min_max['max'] = max(codon_values)

        average_codon_freq_per_aa_dict[aa] = aa_avg_min_max

    return average_codon_freq_per_aa_dict


def calculate_min_max(protein, window_size, cu_table, sliding_size=1):
    """

    :param protein: the protein on which calculate the minmax profile
    :param window_size: the window size used to calculate minmax
    :param cu_table: the cu table of the species in which the protein is expressed
    :param sliding_size: the sliding size to move the window
    :return: a list of value representing the minmax profile
    """
    min_max = []
    average_min_max_for_each_aa = average_min_max_codon_freq_per_aa(cu_table)
    codon_list = protein.get_codon_list()

    # initialize the beginning of the sequence with 0
    for i in range(int(window_size / 2)):
        min_max.append(0)

    # start minimax calculation
    for i in range(0, len(protein.get_codon_list()) - window_size, sliding_size):
        # extract the codons of the window to be analyzed
        window = codon_list[i:i + window_size]

        # get the average usage frequency for each aa and the current usage frequency of the window
        x_avg = sum([average_min_max_for_each_aa[CODON_TO_AA_DICT[codon]]['average'] for codon in window]) / window_size
        x = sum([cu_table.get_codon_usage_table()[codon] for codon in window]) / window_size

        if x > x_avg:
            # if the current is above the average, compute Max value and add it to the list
            x_max = sum([average_min_max_for_each_aa[CODON_TO_AA_DICT[codon]]['max'] for codon in window]) / window_size
            per_max = round((x - x_avg) / (x_max - x_avg) * 100, 2)
            min_max.append(per_max)
        else:
            # else compute Min value and add it to the list
            x_min = sum([average_min_max_for_each_aa[CODON_TO_AA_DICT[codon]]['min'] for codon in window]) / window_size
            per_min = - round(((x_avg - x) / (x_avg - x_min)) * 100, 2)
            min_max.append(per_min)

    # fills in values for codons where window size makes min/max unable to be calculated
    if window_size % 2 == 1:
        for i in range(int(window_size / 2)):
            min_max.append(0)
    else:
        for i in range(int(window_size / 2) - 1):
            min_max.append(0)

    return min_max

import pandas as pd
from scipy.stats import pearsonr


CODON_SIZE = 3

AA_TO_CODON_DICT = {
    'S': ['TCA', 'TCC', 'TCT', 'TCG', 'AGC', 'AGT'],
    'R': ['AGA', 'CGT', 'CGA', 'AGG', 'CGG', 'CGC'],
    'L': ['TTG', 'CTC', 'TTA', 'CTT', 'CTA', 'CTG'],
    'P': ['CCG', 'CCT', 'CCA', 'CCC'],
    'T': ['ACT', 'ACG', 'ACC', 'ACA'],
    'G': ['GGG', 'GGC', 'GGT', 'GGA'],
    'V': ['GTC', 'GTG', 'GTA', 'GTT'],
    'A': ['GCA', 'GCT', 'GCC', 'GCG'],
    '*': ['TGA', 'TAA', 'TAG'],  # Stop codons
    'I': ['ATC', 'ATA', 'ATT'],
    'N': ['AAT', 'AAC'],
    'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'],
    'F': ['TTC', 'TTT'],
    'H': ['CAC', 'CAT'],
    'K': ['AAG', 'AAA'],
    'Y': ['TAT', 'TAC'],
    'C': ['TGC', 'TGT'],
    'Q': ['CAG', 'CAA'],
    'W': ['TGG'],
    'M': ['ATG']
}

CODON_TO_AA_DICT = {'TCA': 'S', 'AAT': 'N', 'TGG': 'W', 'GAT': 'D', 'GAA': 'E', 'TTC': 'F', 'CCG': 'P',
                    'ACT': 'T', 'GGG': 'G', 'ACG': 'T', 'AGA': 'R', 'TTG': 'L', 'GTC': 'V', 'GCA': 'A',
                    'TGA': '*', 'CGT': 'R', 'CAC': 'H', 'CTC': 'L', 'CGA': 'R', 'GCT': 'A', 'ATC': 'I',
                    'ATA': 'I', 'TTT': 'F', 'TAA': '*', 'GTG': 'V', 'GCC': 'A', 'GAG': 'E', 'CAT': 'H',
                    'AAG': 'K', 'AAA': 'K', 'GCG': 'A', 'TCC': 'S', 'GGC': 'G', 'TCT': 'S', 'CCT': 'P',
                    'GTA': 'V', 'AGG': 'R', 'CCA': 'P', 'TAT': 'Y', 'ACC': 'T', 'TCG': 'S', 'ATG': 'M',
                    'TTA': 'L', 'TGC': 'C', 'GTT': 'V', 'CTT': 'L', 'CAG': 'Q', 'CCC': 'P', 'ATT': 'I',
                    'ACA': 'T', 'AAC': 'N', 'GGT': 'G', 'AGC': 'S', 'CGG': 'R', 'TAG': '*', 'CGC': 'R',
                    'AGT': 'S', 'CTA': 'L', 'CAA': 'Q', 'CTG': 'L', 'GGA': 'G', 'TGT': 'C', 'TAC': 'Y',
                    'GAC': 'D'}

def delta_min_max(minmax1, minmax2):
    '''
    :param minmax1: minmax profile of protein 1 (a list of values)
    :param minmax2: minmax profile of protein 2 (a list of values)
    :return: the sum of the  absolute difference between each value of both lists divided by the number of values
    '''
    if len(minmax1) != len(minmax2):
        raise ValueError('Both MinMax profiles must have the same length')

    absolute_differences = []
    for i in range(len(minmax1)):
        absolute_differences.append(abs(minmax1[i]-minmax2[i]))

    return sum(absolute_differences)/len(absolute_differences)

def correlation_min_max(minmax1, minmax2):
    '''
    :param minmax1: minmax profile of protein 1 (a list of values)
    :param minmax2: minmax profile of protein 2 (a list of values)
    :return: the pearson correlation between minmax1 and minmax2 values
    '''
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
    min_max = []
    average_min_max_for_each_aa = average_min_max_codon_freq_per_aa(cu_table)
    codon_list = protein.get_codon_list()
    for i in range(int(window_size / 2)):
        min_max.append(0)

    for i in range(0, len(protein.get_codon_list()) - window_size, sliding_size):
        window = codon_list[i:i + window_size]

        x_avg = sum([average_min_max_for_each_aa[CODON_TO_AA_DICT[codon]]['average'] for codon in window]) / window_size
        x = sum([cu_table.get_codon_usage_table()[codon] for codon in window]) / window_size

        if x > x_avg:
            x_max = sum([average_min_max_for_each_aa[CODON_TO_AA_DICT[codon]]['max'] for codon in window]) / window_size
            per_max = round((x - x_avg) / (x_max - x_avg) * 100, 2)
            min_max.append(per_max)
        else:
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


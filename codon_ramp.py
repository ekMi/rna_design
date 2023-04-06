import random
from utils.codon_aa_dict import CODON_TO_AA_DICT, AA_TO_CODON_DICT
from codon_usage.codon_usage import CodonUsage
from utils.protein import Protein


def replace_codons_with_rare(sequence: str, cu_table: CodonUsage, threshold=0.4, length_to_change=16):
    """
    :param cu_table: the codon usage of the species to define the codon frequencies
    :param sequence: the inputu protein on which to reduce the codon usage
    :param threshold: the threshold to define a codon as rare
    :param length_to_change: the number of codon to change
    :return: the modified protein with codon freq bellow threshold but above 0.1
    """

    frequencies = cu_table.get_relative_codon_freq_table()
    seq = sequence[:length_to_change*3]

    # for each codon
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        # if the current codon is already below the threshold, keep it
        if 0.1 < frequencies[codon] < threshold:
            new_codon = codon
        else:
            # Get the codon possibilities that have a usage frequency above 0.1, but below the threshold
            aa = CODON_TO_AA_DICT[codon]
            possible_codons = AA_TO_CODON_DICT[aa]
            rare_codons = [c for c in possible_codons if 0.1 < frequencies[c] < threshold]
            if rare_codons:
                new_codon = random.choice(rare_codons)
            # if the constraint returns an empty list, get the codon with the lowest frequency but still over 0.1
            else:
                alt = [c for c in AA_TO_CODON_DICT[aa] if 0.1 < frequencies[c]]
                new_codon = min(alt, key=lambda c: frequencies[c])

        seq = seq[:i] + new_codon + seq[i + 3:]

    return seq



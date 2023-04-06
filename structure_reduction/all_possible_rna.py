from itertools import product
from utils.codon_aa_dict import CODON_TO_AA_DICT, AA_TO_CODON_DICT
from structure_reduction.struct_reduction import StructReductionStrategy, evaluate_structure
import random


def generate_all_possible_mrnas(sd, spacer, cds, cu_table, threshold):
    """
    This function generates all possible mRNA alternatives that code for the same protein.

    :param sd: str, the Shine-Dalgarno sequence
    :param spacer: str, the spacer sequence between the SD and CDS
    :param cds: str, the coding region of the RNA
    :param cu_table: codon usage table, an object that holds the relative frequency of each codon
    :param threshold: float, the maximum allowed frequency for a codon to be considered rare
    :return: list, a list with all possible mRNA alternatives coding for the same protein
    """

    # Create a list of ACGT for the number of nucleotides of the spacer
    possible_nct = [['A', 'C', 'G', 'T'] for i in range(len(spacer))]

    # Create all the combinations for the spacer
    spacers = list(product(*possible_nct))

    # Create all possible combinations for the given cds size with the frequency threshold constraint
    cds_to_codon = [cds[i:i + 3] for i in range(0, len(cds), 3)]  # split CDS into codons
    amino_acid_sequence = [CODON_TO_AA_DICT[codon] for codon in cds_to_codon]  # get amino acid sequence from codons
    frequencies = cu_table.get_relative_codon_freq_table()  # get codon usage frequencies
    possible_codons = []
    for i in range(10):
        amino_acid = amino_acid_sequence[i]
        rare_codons = [c for c in AA_TO_CODON_DICT[amino_acid] if 0.1 < frequencies[c] < threshold]  # rare codons
        if rare_codons:
            allowed_codons = rare_codons  # if there are rare codons, only allow those
        else:
            # if no codon found between 0.1 and threshold, take the lowest (still over 0.1)
            alt = [c for c in AA_TO_CODON_DICT[amino_acid] if 0.1 < frequencies[c]]  # alternative codons
            allowed_codons = [min(alt, key=lambda c: frequencies[c])]  # select the least frequent alternative codon
        possible_codons.append(allowed_codons)
    all_possible_cds = list(product(*possible_codons))  # create all possible combinations of codons

    # Create all combinations of spacers * possible CDS
    sequences = []
    for combination in all_possible_cds:
        modified_codons = ''.join(combination)  # join codons into a single string
        modified_cds = modified_codons + cds[10 * 3:]  # modify CDS by replacing first 10 codons with selected ones
        for spacer in spacers:
            modified_spacer = ''.join(spacer)
            sequence = sd + modified_spacer + modified_cds  # create final sequence
            sequences.append(sequence)

    random.shuffle(sequences)
    return sequences


class TestAllRnasStrategy(StructReductionStrategy):

    def run(self, sequence, sd, spacer, cu_table, threshold=0.4):
        print('Starting all possible rnas optimization')
        print('---------------------------------------')
        all_possible_mrnas = generate_all_possible_mrnas(sd, spacer, sequence, cu_table, threshold)
        print(f"Number of possible sequences = {len(all_possible_mrnas)}")

        best_seq = ''
        best_score = len(sd + spacer + sequence)
        n_best_score = 0
        iteration_number_best = 0
        for i, mrna in enumerate(all_possible_mrnas):
            new_score, new_structure = evaluate_structure(mrna, len(sd + spacer) + 30)
            if new_score < best_score:
                best_score = new_score
                best_seq = mrna
                n_best_score = 0
                iteration_number_best = i
            if new_score == best_score:
                n_best_score += 1
            if i % 500 == 0:
                print(f'{i} solutions tested')

        best_spacer = best_seq[len(sd):len(sd + spacer)]
        best_sequence = best_seq[len(sd + spacer):]
        print(f'Number of sequences with best score: {n_best_score}')
        return best_spacer, best_sequence, best_score, iteration_number_best

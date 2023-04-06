from structure_reduction.struct_reduction import StructReductionStrategy, evaluate_structure
from utils.codon_aa_dict import CODON_TO_AA_DICT, AA_TO_CODON_DICT
from itertools import product
from codon_usage.codon_usage import CodonUsage
import random
import math


def count_all_possible_mrnas(spacer, cds, cu_table, threshold):
    """
    This function counts all possible mRNA alternatives that code for the same protein.

    :param spacer: str, the spacer sequence between the SD and CDS
    :param cds: str, the coding region of the RNA
    :param cu_table: codon usage table, an object that holds the relative frequency of each codon
    :param threshold: float, the maximum allowed frequency for a codon to be considered rare
    :return: int, the number of possible mRNA alternatives coding for the same protein
    """

    # Count all possible combinations for the given cds size with the frequency threshold constraint
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
    num_possible_cds = len(list(product(*possible_codons)))  # count all possible combinations of codons

    # Count all possible combinations of spacers * possible CDS
    num_possible_spacers = 4**len(spacer)  # count all possible combinations of spacers
    num_possible_mrnas = num_possible_spacers * num_possible_cds

    return num_possible_mrnas


def get_alternative_codon(aa: str, cu_table: CodonUsage, current_codon: str, threshold):
    frequencies = cu_table.get_relative_codon_freq_table()
    possible_codons = AA_TO_CODON_DICT[aa]
    rare_codons = [c for c in possible_codons if 0.1 < frequencies[c] < threshold]

    if len(rare_codons) > 1 and current_codon in rare_codons:
        rare_codons.remove(current_codon)

    if rare_codons:
        new_codon = random.choice(rare_codons)
    # if the constraint returns an empty list, get the codon with the lowest frequency but still over 0.1
    else:
        alt = [c for c in AA_TO_CODON_DICT[aa] if 0.1 < frequencies[c]]
        new_codon = min(alt, key=lambda c: frequencies[c])

    return new_codon


def mutate_sequence_random(spacer, sequence, current_temperature, init_temperature, cu_table, threshold):
    codons = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
    # calculate the number of mutations to perform
    n_mutations = int(10 * (current_temperature / init_temperature))

    # randomly mutate the spacer
    spacer_list = list(spacer)
    id_changed = []
    for i in range(min(n_mutations, len(spacer_list))):
        idx = random.randint(0, len(spacer_list) - 1)
        while idx in id_changed:
            idx = random.randint(0, len(spacer_list) - 1)
        id_changed.append(idx)
        spacer_list[idx] = random.choice(['A', 'C', 'G', 'T'])
    new_spacer = ''.join(spacer_list)

    # mutate the first codons
    id_changed = []
    for i in range(n_mutations):
        idx = random.randint(0, 9)
        # make sure we don't chose an id already changed for this mutation
        while idx in id_changed:
            idx = random.randint(0, 9)
        id_changed.append(idx)

        aa = CODON_TO_AA_DICT[codons[idx]]

        new_codon = get_alternative_codon(aa, cu_table, codons[idx], threshold)
        codons[idx] = new_codon

    # rejoin the sequence and return it
    return new_spacer, ''.join(codons)


def is_accepted_new_solution(current, new, temperature):
    if new < current:
        # if the new score is better, the solution is accepted
        return True
    else:
        # else use a Boltzmann distribution to accept or not the solution
        # The higher the temperature, the higher the probability to accept a less optimal solution
        return random.random() < math.exp((new - current) / temperature)


def update_temperature(iteration_number, max_iterations, initial_temperature, current_temperature):
    # cooling_rate = current_temperature / initial_temperature
    cooling_rate = 0.9
    iteration_ratio = iteration_number / max_iterations
    new_temperature = initial_temperature * (1 - cooling_rate * iteration_ratio)
    return new_temperature


class SimulatedAnnealing(StructReductionStrategy):
    def __init__(self, max_iteration=1000, max_time=3600, initial_temperature=100.0, cooling_rate=0.9):
        self._max_time = max_time
        self._initial_temperature = initial_temperature
        self._cooling_rate = cooling_rate
        self._max_iteration = max_iteration

    def run(self, sequence, sd, spacer, cu_table, threshold=0.4):
        print('Starting simulated annealing optimization')
        print('_________________________________________')
        n_nucleotides_to_evaluate = len(sd + spacer + sequence[:30])

        current_seq = sequence[:]
        current_spacer = spacer[:]
        current_seq_score = evaluate_structure(sd + current_spacer + current_seq, n_nucleotides_to_evaluate)[0]
        print(f'Init score is {current_seq_score}')

        best_seq = current_seq
        best_spacer = current_spacer
        best_score = current_seq_score
        current_temperature = self._initial_temperature
        iteration_number_best = 0
        for i in range(self._max_iteration):
            # Generate new sequence
            new_spacer, new_seq = mutate_sequence_random(current_spacer, current_seq, current_temperature,
                                                         self._initial_temperature, cu_table, threshold)

            # Calculate its score
            new_score = evaluate_structure(sd + new_spacer + new_seq, n_nucleotides_to_evaluate)[0]

            # Check if new solution is accepted or not
            if is_accepted_new_solution(current_seq_score, new_score, current_temperature):
                current_seq = new_seq
                current_spacer = new_spacer
                current_seq_score = new_score

            # Update of the best solution found so far
            if new_score < best_score:
                best_score = new_score
                best_seq = new_seq
                best_spacer = new_spacer
                iteration_number_best = i

            # Update temperature value
            current_temperature = update_temperature(i, self._max_iteration, self._initial_temperature,
                                                     current_temperature)

            # if a solution is found (no structure for the start of the sequence, stop iteration
            if best_score == 0:
                break

        return best_spacer, best_seq, best_score, iteration_number_best

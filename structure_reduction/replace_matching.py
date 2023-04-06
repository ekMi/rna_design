from utils.codon_aa_dict import CODON_TO_AA_DICT, AA_TO_CODON_DICT
from codon_usage.codon_usage import CodonUsage
import random
from structure_reduction.struct_reduction import StructReductionStrategy, evaluate_structure

SCORES = {
    'GC': 3.12,
    'CG': 3.12,
    'AT': 1,
    'TA': 1,
    'GT': 1,
    'TG': 1
}

MATCHING = {
    "A": ["T"],
    "T": ["A", "G"],
    "G": ["C", "T"],
    "C": ["G"]
}


def get_closing_parenthesis_position(seq, start):
    """
    :param seq: a dot parenthesis sequence
    :param start: the position to start the exploration (the opening parenthesis + 1)
    :return: the position of the corresponding closing parenthesis
    """
    opened_parenthesis = 1
    for i in range(start + 1, len(seq)):
        if seq[i] == '(':
            opened_parenthesis += 1
        elif seq[i] == ')':
            opened_parenthesis -= 1

        if opened_parenthesis == 0:
            return i


def matching_score(codon1, codon2):
    """
    Calculates the matching score between two codons. The matching score is the sum of the scores of matching
    letters at each position in the codons.

    :param codon1: A 3-letter string representing a codon
    :param codon2: A 3-letter string representing a codon
    :return: An integer representing the matching score between the two codons
    """
    score = 0
    # Loop over each position in the codons
    for i in range(3):
        # Create a string of the two letters to match at this position
        matching_letters = f'{codon1[i]}{codon2[i]}'
        # Add the score for this matching pair to the total score
        score += SCORES.get(matching_letters, 0)
    # Return the total score for the two codons
    return score


def get_alternatives(codon: str, matching_codon: str, matching_replacement_allowed: bool, cu_table: CodonUsage,
                     first_threshold=0.4, second_threshold=1.0):
    """
    This function returns a tuple with the best alternative codon and matching codon, given a starting codon, a matching
    codon, and a codon usage table. The function optimizes for codons that are rare but still have a high enough frequency
    in the codon usage table. If matching_replacement_allowed is set to True, the matching codon can also be changed to
    another rare but frequently used codon.

    :param matching_replacement_allowed: boolean value that states if the matching codon can be changed
    :param codon: the codon originating from the start of the sequence
    :param matching_codon: the codon for which there is a matching
    :param cu_table: CodonUsage object containing codon usage information
    :param threshold: maximum frequency of a codon to be considered "rare"
    :return: a tuple: (codon, matching codon) as a better alternative to reduce the structure
    """

    # Get relative codon frequency table from codon usage object
    frequencies = cu_table.get_relative_codon_freq_table()

    # Get rare codons that are alternatives to the starting codon
    codon_alternatives = [c for c in AA_TO_CODON_DICT[CODON_TO_AA_DICT[codon]] if
                          0.1 < frequencies[c] < first_threshold]

    # If there are rare codons, only allow those
    if codon_alternatives:
        allowed_codons = codon_alternatives
    # If no rare codon is found, take the lowest frequency codon that is still over 0.1
    else:
        alt = [c for c in AA_TO_CODON_DICT[CODON_TO_AA_DICT[codon]] if 0.1 < frequencies[c]]
        allowed_codons = [min(alt, key=lambda c: frequencies[c])]

    # Get alternatives for the matching codon
    if matching_replacement_allowed:
        matching_codon_alternatives = [c for c in AA_TO_CODON_DICT[CODON_TO_AA_DICT[matching_codon]] if
                                       0.1 < frequencies[c] < second_threshold]
        if matching_codon_alternatives:
            allowed_matching_codons = matching_codon_alternatives
        else:
            alt = [c for c in AA_TO_CODON_DICT[CODON_TO_AA_DICT[matching_codon]] if 0.1 < frequencies[c]]
            allowed_matching_codons = [min(alt, key=lambda c: frequencies[c])]
    else:
        allowed_matching_codons = [matching_codon]

    best_score = 50
    best_codon_alt = ""
    best_matching_alt = ""
    for codon_alt in allowed_codons:
        for matching_alt in allowed_matching_codons:
            score = matching_score(codon_alt, matching_alt)
            if score < best_score:
                best_score = score
                best_codon_alt = codon_alt
                best_matching_alt = matching_alt
    return best_codon_alt, best_matching_alt


def replace_matching_sd(nucleotide, matching_position, sd, spacer, seq):
    """
    Replace a nucleotide in the spacer sequence or a codon in the coding sequence
    that do not match the given nucleotide with a random alternative.

    :param:
    - nucleotide: the nucleotide to match
    - matching_position: the position of the matching nucleotide in the sequence
    - sd: the Shine-Dalgarno sequence
    - spacer: the spacer sequence
    - seq: the coding sequence

    :return:
    - a tuple of the updated spacer and coding sequences
    """
    # if matching is in the spacer
    if len(sd) < matching_position < len(sd + spacer):
        # replace by a single random nuc not matching with nuc
        alternatives = [alt for alt in ["A", "T", "G", "C"] if alt not in MATCHING[nucleotide]]
        alt = random.choice(alternatives)
        pos = matching_position - len(sd)
        new_spacer = spacer[:pos] + alt + spacer[pos:]
        return new_spacer, seq
    # else it is in the coding region
    else:
        pos = matching_position - len(sd + spacer)
        match_first_pos = (pos // 3) * 3
        matching_codon = seq[match_first_pos: match_first_pos + 3]
        codon_alternatives = [codon for codon in AA_TO_CODON_DICT[CODON_TO_AA_DICT[matching_codon]]
                              if codon not in MATCHING[nucleotide]]
        codon_alternative = random.choice(codon_alternatives)
        new_seq = seq[:match_first_pos] + codon_alternative + seq[match_first_pos + 3:]
        return spacer, new_seq


def replace_matching_spacer(spacer: str, nucleotide_position: int, matching_nucleotide: str):
    alternatives = [alt for alt in ["A", "T", "G", "C"] if alt not in MATCHING[matching_nucleotide]]
    alt = random.choice(alternatives)
    new_spacer = spacer[:nucleotide_position] + alt + spacer[nucleotide_position + 1:]
    return new_spacer


def replace_matching_codon(nucleotide_position: int, matching_position: int, sequence: str,
                           n_nucleotides_to_evaluate: int, cu_table: CodonUsage,
                           replacement_outside_first_allowed: bool, threshold: float):
    first_codon_pos = (nucleotide_position // 3) * 3
    match_codon_pos = (matching_position // 3) * 3
    first_codon = sequence[first_codon_pos: first_codon_pos + 3]
    matching_codon = sequence[match_codon_pos: match_codon_pos + 3]

    matching_replacement_allowed = matching_position < n_nucleotides_to_evaluate or replacement_outside_first_allowed
    if match_codon_pos > n_nucleotides_to_evaluate:
        second_threshold = 1
    else:
        second_threshold = threshold
    first_alternative, matching_alternative = get_alternatives(first_codon, matching_codon,
                                                               matching_replacement_allowed, cu_table,
                                                               second_threshold=second_threshold)
    new_seq = sequence[:first_codon_pos] + first_alternative + sequence[
                                                               first_codon_pos + 3: match_codon_pos] + matching_alternative + sequence[
                                                                                                                              match_codon_pos + 3:]
    return new_seq


def mutate_sequence_avoid_matching(sd: str, spacer: str, seq: str, struct: str, n_nucleotides_to_evaluate: int,
                                   cu_table: CodonUsage, replacement_outside_first_allowed: bool, threshold: float):
    # check for opening parenthesis in the beginning
    new_sequence = seq[:]
    new_spacer = spacer[:]
    i = 0
    while i < n_nucleotides_to_evaluate:
        if struct[i] == '(':
            match_position = get_closing_parenthesis_position(struct, i)
            # check if this position is in the SD, the spacer or first codons
            if i < len(sd):
                # this is in the sd, replace only the matching codon in the first 10 only if replacement outside not
                # allowed
                if match_position < n_nucleotides_to_evaluate or replacement_outside_first_allowed:
                    new_spacer, new_sequence = replace_matching_sd(sd[i], match_position, sd, new_spacer, new_sequence)
                i += 1
            elif len(sd) <= i < len(sd + new_spacer):
                # this is in the spacer, replace by an alternative nucleotide not pairing
                new_spacer = replace_matching_spacer(new_spacer, i - len(sd),
                                                     new_sequence[match_position - len(sd + new_spacer)])
                i += 1
            else:
                # this is in the coding region, replace with codons
                new_sequence = replace_matching_codon(i - len(sd + new_spacer), match_position - len(sd + new_spacer),
                                                      new_sequence, n_nucleotides_to_evaluate, cu_table,
                                                      replacement_outside_first_allowed, threshold)
                # Then jump i to next codon
                i = (i // 3) * 3 + 3
        i += 1
    return new_spacer, new_sequence


class ReplaceMatchingCodonsStrategy(StructReductionStrategy):
    def __init__(self, replacement_outside_first_allowed=False):
        self._replacement_outside_first_allowed = replacement_outside_first_allowed

    def run(self, sequence, sd, spacer, cu_table, threshold=0.4):
        if self._replacement_outside_first_allowed:
            message = "with"
        else:
            message = " without"

        print(f'Starting replace matching optimization replacement {message} outside first codons allowed')
        print('__________________________________________________________________________________________')

        rbs = sd + spacer
        n_nucleotides_to_evaluate = len(rbs) + 30

        current_seq = sequence[:]
        current_spacer = spacer[:]
        current_seq_score, current_seq_struct = evaluate_structure(sd + current_spacer + current_seq,
                                                                   n_nucleotides_to_evaluate)

        print(f'Init score is {current_seq_score}')

        # Loop until the score is above 0 or if the score did not evolves after 100 iterations
        flag = 0
        iteration_number = 0
        best_seq = current_seq
        best_spacer = current_spacer
        best_seq_struct = current_seq_struct
        best_score = current_seq_score
        iteration_number_best = 0
        while current_seq_score > 0 and flag < 100:
            # mutate sequence for element involved in a structure
            new_spacer, new_sequence = mutate_sequence_avoid_matching(sd, best_spacer, best_seq, best_seq_struct,
                                                                      n_nucleotides_to_evaluate, cu_table,
                                                                      self._replacement_outside_first_allowed, threshold)
            # calculate score
            new_score, new_struct = evaluate_structure(sd + new_spacer + new_sequence, n_nucleotides_to_evaluate)
            # if best score, keep the mutated sequence for next iteration
            if new_score < best_score:
                best_score = new_score
                best_seq = new_sequence
                best_seq_struct = new_struct
                best_spacer = new_spacer
                flag = 0
                iteration_number_best = iteration_number
            else:
                flag += 1

            iteration_number += 1

        return best_spacer, best_seq, best_score, iteration_number_best

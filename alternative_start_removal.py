import re
import random
from utils.codon_aa_dict import CODON_TO_AA_DICT, AA_TO_CODON_DICT, CODON_SIZE
from utils.import_proteins_file import import_proteins

def find_alternative_starts(seq: str, sd_alternatives: list = ['AGGAGG']):
    """
    Find all alternative start codon indices in a given nucleotide sequence.

    Args:
        seq (str): nucleotide sequence in string format
        sd_alternatives (list): a list of Shine-Dalgarno sequence alternatives to search for.
                                Default is set to ['AGGAGG']

    Returns:
        list: a list of all alternative start codon indices found in the sequence
    """
    alternative_start_indices = []
    for alt in sd_alternatives:
        pattern = f'{alt}[A,T,G,C]{{4,14}}ATG'  # create a regular expression pattern to find alternative start sites
        matches = re.finditer(pattern, seq)     # find all matches for the pattern in the sequence
        # collect all alternative start indices. Skip the first one as it is the real start
        indices = [match.start() for match in matches if match != 0]
        alternative_start_indices += indices
    return alternative_start_indices


def replace_codons(seq, indices):
    """
    Replace the codons at the given indices with a random alternative codon.

    Args:
        seq (str): nucleotide sequence in string format
        indices (list): a list of indices where alternative start codons are to be replaced

    Returns:
        str: a new nucleotide sequence with alternative start codons replaced
    """
    seq_copy = seq[:]  # create a copy of the input sequence to avoid modifying the original sequence
    for idx in indices:
        codon_start_pos = (idx // CODON_SIZE) * CODON_SIZE  # find the starting position of the codon
        codon_to_change = seq_copy[codon_start_pos : codon_start_pos + CODON_SIZE]  # extract the codon to be changed
        codon_alternatives = AA_TO_CODON_DICT[CODON_TO_AA_DICT[codon_to_change]]  # find all alternative codons for the amino acid
        codon_alternatives.remove(codon_to_change)  # remove the original codon from the list of alternatives
        replace_codon = random.choices(codon_alternatives, k=1)[0]  # randomly select an alternative codon
        seq_copy = seq_copy[:codon_start_pos] + replace_codon + seq_copy[codon_start_pos+3:]  # replace the original codon with the new one
    return seq_copy  # return the modified sequence


if __name__ == '__main__':
    proteins = import_proteins(file_path='./test/test_files/test_proteins.txt')
    rpoD = proteins[2]
    test_seq = rpoD.get_dna_seq()+'AGGAGG'+'AAAAAAA'+'ATG'
    print(test_seq[-20:])
    alt_starts = find_alternative_starts(test_seq)
    new_seq = replace_codons(test_seq, alt_starts)
    print(new_seq[-20:])

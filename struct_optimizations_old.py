from utils.protein import Protein
import os

N_NUCLEOTIDES_STRUCT_REDUCTION = 30

SCORES = {
    'GC': 3.12,
    'CG': 3.12,
    'AU': 1,
    'UA': 1,
    'GU': 1,
    'UG': 1
}

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


def get_closing_parenthesis_position(seq, start):
    """
    :param seq: a dot parenthesis sequence
    :param start: the position to start the exploration (the opening parenthesis + 1)
    :return: the position of the corresponding closing parenthesis
    """
    opened_parenthesis = 1
    for i in range(start, len(seq)):
        if seq[i] == '(':
            opened_parenthesis += 1
        elif seq[i] == ')':
            opened_parenthesis -= 1

        if opened_parenthesis == 0:
            return i


def matching_score(codon1, codon2):
    score = 0
    for i in range(3):
        matching_letters = f'{codon1[i]}{codon2[i]}'
        score = score + SCORES.get(matching_letters, 0)
    return score


def get_alternatives(codon, matching_codon):
    """
    :param codon: the codon originating from the start of the sequence
    :param matching_codon: the codon for which there is a matching
    :return: a tuple: (codon, matching codon) as a better alternative to reduce the structure
    """
    codon_alternatives = AA_TO_CODON_DICT[CODON_TO_AA_DICT[codon]]
    matching_codon_alternatives = AA_TO_CODON_DICT[CODON_TO_AA_DICT[matching_codon]]

    best_score = 10
    best_codon_alt = ""
    best_matching_alt = ""
    for codon_alt in codon_alternatives:
        for matching_alt in matching_codon_alternatives:
            score = matching_score(codon_alt, matching_alt)
            if score < best_score:
                best_score = score
                best_codon_alt = codon_alt
                best_matching_alt = matching_alt
    return best_codon_alt, best_matching_alt


def struct_reduction(seq, sd, spacer):
    """
    :param seq: the rna sequence for which the structure needs to be reduced for the firsts 30 nucleotides
    :return: the rna sequence for which the structure has been reduced for the 30 first nucleotides
    """

    # initialize with the given structure and add  SD and spacer
    rbs = sd + spacer
    seq = rbs + seq

    command = f'echo {seq} | RNAfold -C'
    rna_fold_out = os.popen(command).read()
    seq_struct = rna_fold_out.split('\n')[1].split(' ')[0]
    print('origin')
    print(seq_struct)
    alternative_struct = []
    alternative_seqs = []
    flag = 0
    while flag < 10:
        for j in range(len(rbs), N_NUCLEOTIDES_STRUCT_REDUCTION+len(rbs)):
            if seq_struct[j] == '(':
                # if a structure is found, get the matching codons
                codon_first_nuc_position = (j // 3) * 3
                first_codon = seq[codon_first_nuc_position: codon_first_nuc_position + 3]

                closing_parenthesis_position = get_closing_parenthesis_position(seq_struct, j + 1)
                matching_codon_first_nuc_position = (closing_parenthesis_position // 3) * 3
                matching_codon = seq[matching_codon_first_nuc_position:matching_codon_first_nuc_position + 3]

                codon_alternative, matching_alternative = get_alternatives(first_codon, matching_codon)

                # replace the codons at their respective positions in the seq
                alternative_seq = seq[:codon_first_nuc_position]
                alternative_seq += codon_alternative
                alternative_seq += seq[codon_first_nuc_position + 3:matching_codon_first_nuc_position]
                alternative_seq += matching_alternative
                alternative_seq += seq[matching_codon_first_nuc_position + 3:]

                seq = alternative_seq

                # test if it improves the situation
        command = f'echo {seq} | RNAfold -C'
        rna_fold_out = os.popen(command).read()
        optimized_struc = rna_fold_out.split('\n')[1].split(' ')[0]
        alternative_struct.append(optimized_struc[:N_NUCLEOTIDES_STRUCT_REDUCTION+len(rbs)].count('('))
        alternative_seqs.append(seq)

        flag += 1

    # find the best alternative struct (the one with less '(')
    best = N_NUCLEOTIDES_STRUCT_REDUCTION+len(rbs)
    for i, alt in enumerate(alternative_struct):
        if alt < best:
            print(alt)
            best = alt
            seq = alternative_seqs[i]

    return seq


if __name__ == '__main__':

    protein_file = open('test/test_files/test_proteins.txt', 'r')

    lines = protein_file.readlines()

    proteins = []
    for count, line in enumerate(lines):
        if ":" in line:
            if count != 0:
                proteins.append(Protein(protein_name, sequence, 'dna'))
            sequence = ""
            protein_name = line.strip()
        else:
            sequence += line.strip()
    proteins.append(Protein(protein_name, sequence, 'dna'))

    rpoD = proteins[2]

    sd = 'AGGAGG'
    spacer = ''
    improved_rpoD = struct_reduction(rpoD.get_dna_seq(), sd, spacer)
    print(sd+spacer+rpoD.get_dna_seq()[:99])
    print(improved_rpoD)

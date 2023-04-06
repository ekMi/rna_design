from utils.protein import Protein
from codon_usage.cai import calculate_cai
from codon_usage.codon_usage import CodonUsage
import random
from utils.codon_aa_dict import AA_TO_CODON_DICT


def get_random_codons_pool_for_each_aa(codon_freq_table, count_of_each_aa):
    codon_pools_per_aa = dict()
    # for each aa, get the number of requested codons
    for aa, count in count_of_each_aa.items():
        codon_pools_per_aa[aa] = []
        # for each codon related to this aa, creates the number of codons related to the frequency of use from cu table
        for codon in AA_TO_CODON_DICT[aa]:
            codon_freq = codon_freq_table[codon]
            number_of_codons = round(codon_freq * count)
            for i in range(0, number_of_codons):
                codon_pools_per_aa[aa].append(codon)
        # due to rounding errors, it is possible that 1 codon could be missed
        # in such case, add the most frequent codon to the list
        while len(codon_pools_per_aa[aa]) < count:
            codons = AA_TO_CODON_DICT[aa]
            codons_freq = []
            for codon in codons:
                codons_freq.append(codon_freq_table[codon])
            most_frequent_codon = [x for _, x in sorted(zip(codons_freq, codons), reverse=True)][0]
            codon_pools_per_aa[aa].append(most_frequent_codon)
        # shuffle the codons of the aa
        random.shuffle(codon_pools_per_aa[aa])

    return codon_pools_per_aa


def randomize(protein: Protein, cu_table: CodonUsage):
    """
    :param protein: the protein on wich to randomize the codon usage keeping it's aa sequence
    :param cu_table: the cu table of the species in which the protein is expressed
    :return: a protein in which the codon choices have been randomized
                taking into account the codon usage of the species
    """
    # get a pool of codons for each amino acid
    pool_of_codons = get_random_codons_pool_for_each_aa(cu_table.get_relative_codon_freq_table(),
                                                        protein.amino_acid_counts)
    randomized_sequence = ''
    # for each amino acid of the sequence, choose one codon in the list of codons
    for aa in protein.get_aa_seq():
        randomized_sequence += pool_of_codons[aa].pop(0)

    # recreate the protein object to return
    randomized_prot = Protein(protein.get_name(), randomized_sequence, 'dna')
    return randomized_prot


if __name__ == '__main__':
    cu = CodonUsage("../test/test_files/CUB_HiveCut_Ecoli_K12.txt", "cub")
    rpsu_dna_sequence = 'ATGCCGGTAATTAAAGTACGTGAAAACGAGCCGTTCGACGTAGCTCTGCGTCGCTTCAAGCGTTCCTGCGAAAAAGCAGGTGTTCTGGCGGAAGTTCGTCGTCGTGAGTTCTATGAAAAACCGACTACCGAACGTAAGCGCGCTAAAGCTTCTGCAGTGAAACGTCACGCGAAGAAACTGGCTCGCGAAAACGCACGCCGCACTCGTCTGTACTAA'
    rpsu = Protein('rpsu', rpsu_dna_sequence, 'dna')
    rpsu_randomized = randomize(rpsu, cu)

    cu_ser = CodonUsage("../test/test_files/ScerCUB.txt", "cub")
    rpsu_ser_rando = randomize(rpsu, cu_ser)
    print(calculate_cai(rpsu, cu_ser))
    print(calculate_cai(rpsu_randomized, cu_ser))

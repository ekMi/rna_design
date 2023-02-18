from protein import Protein
from cai import calculate_cai
from codon_usage import CodonUsage
import random

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

def get_random_codons_pool_for_each_aa(codon_freq_table, count_of_each_aa):
    codon_pools_per_aa = dict()
    for aa, count in count_of_each_aa.items():
        codon_pools_per_aa[aa] = []
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
        random.shuffle(codon_pools_per_aa[aa])

    return codon_pools_per_aa


def randomize(prot, cu):
    pool_of_codons = get_random_codons_pool_for_each_aa(cu.get_relative_codon_freq_table(), prot.amino_acid_counts)
    randomized_sequence = ''
    for aa in prot.get_aa_seq():
        randomized_sequence += pool_of_codons[aa].pop(0)

    randomized_prot = Protein(prot.get_name(), randomized_sequence, 'dna')
    return randomized_prot


if __name__ == '__main__':
    cu = CodonUsage("test_files/CUB_HiveCut_Ecoli_K12.txt", "cub")
    rpsu_dna_sequence = 'ATGCCGGTAATTAAAGTACGTGAAAACGAGCCGTTCGACGTAGCTCTGCGTCGCTTCAAGCGTTCCTGCGAAAAAGCAGGTGTTCTGGCGGAAGTTCGTCGTCGTGAGTTCTATGAAAAACCGACTACCGAACGTAAGCGCGCTAAAGCTTCTGCAGTGAAACGTCACGCGAAGAAACTGGCTCGCGAAAACGCACGCCGCACTCGTCTGTACTAA'
    rpsu = Protein('rpsu', rpsu_dna_sequence, 'dna')
    rpsu_randomized = randomize(rpsu, cu)


    cu_ser = CodonUsage("test_files/ScerCUB.txt", "cub")
    rpsu_ser_rando = randomize(rpsu, cu_ser)
    print(calculate_cai(rpsu, cu_ser))
    print(calculate_cai(rpsu_randomized, cu_ser))


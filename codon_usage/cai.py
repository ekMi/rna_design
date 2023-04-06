import math
from codon_usage.codon_usage import CodonUsage
from utils.protein import Protein
from utils.codon_aa_dict import CODON_SIZE


def remove_single_codons(seq):
    """
    Remove the ATG and TGG codons since they are the only codon for their AA and therefore do not impact
    the CAI computation
    :param seq: a DNA coding sequence for one specific protein
    :return: the same sequence without{'ATG', 'TGG'}
    """
    cleaned_sequence = ''
    for i in range(0, len(seq), CODON_SIZE):
        if seq[i: i + CODON_SIZE] != 'ATG' and seq[i: i + CODON_SIZE] != 'TGG':
            cleaned_sequence += seq[i: i + CODON_SIZE]
    return cleaned_sequence


def calculate_cai(prot: Protein, cu_table: CodonUsage):
    """
    Calculate the codon adaptation index of a given protein based on a given codon usage table
    :param prot: a protein object
    :param cu_table: the codon usage table (dictionary)
    :return: the codon adaptation index
    """

    rscu_table = cu_table.get_rscu_table()
    # clean the protein from unnecessary codons
    cleaned_prot = Protein(prot.get_name(), remove_single_codons(prot.get_dna_seq()), "dna")

    cleaned_prot_dna_seq = cleaned_prot.get_dna_seq()
    # compute the cai of the current sequence
    freq_obs = []
    for i in range(0, cleaned_prot.nucleotides_count, CODON_SIZE):
        if rscu_table[cleaned_prot_dna_seq[i: i + CODON_SIZE]] > 0:
            freq_obs.append(rscu_table[cleaned_prot_dna_seq[i: i + CODON_SIZE]])
        else:
            freq_obs.append(0.5)
    cai_obs = math.pow(math.prod(freq_obs), 1 / cleaned_prot.total_amino_acids_count)

    # compute the cai based on the assumption of the use of the preferred codon for each AA
    freq_max = []
    cleaned_prot_aa_seq = cleaned_prot.get_aa_seq()
    preferred_codon_rscu_for_each_aa = cu_table.get_preferred_codon_rscu_for_each_aa()
    for i in range(0, cleaned_prot.total_amino_acids_count, 1):
        freq_max.append(preferred_codon_rscu_for_each_aa[cleaned_prot_aa_seq[i]])
    cai_max = math.pow(math.prod(freq_max), 1 / cleaned_prot.total_amino_acids_count)

    # return the ratio of the observed CAI/the maximum CAI
    return round(cai_obs / cai_max, 3)


if __name__ == '__main__':
    cu_table = CodonUsage("../test/test_files/VeryHighLevelExpressionEColi.txt")
    rpsu_prot = Protein('rspu',
                        'ATGCCGGTAATTAAAGTACGTGAAAACGAGCCGTTCGACGTAGCTCTGCGTCGCTTCAAGCGTTCCTGCGAAAAAGCAGGTGTTCTGGCGGAAGTTCGTCGTCGTGAGTTCTATGAAAAACCGACTACCGAACGTAAGCGCGCTAAAGCTTCTGCAGTGAAACGTCACGCGAAGAAACTGGCTCGCGAAAACGCACGCCGCACTCGTCTGTACTAA',
                        'dna')
    print(calculate_cai(rpsu_prot, cu_table))

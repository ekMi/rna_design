import math
from codon_usage import CodonUsage
from protein import Protein

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


def remove_single_codons(seq):
    """
    Remove the ATG and TGG codons since they are the only codon for their AA and therefore do not impact
    the CAI computation
    :param seq: a DNA coding sequence for one specific protein
    :return: the same sequence \{'ATG', 'TGG'}
    """
    cleaned_sequence = ''
    for i in range(0, len(seq), CODON_SIZE):
        if seq[i: i + CODON_SIZE] != 'ATG' and seq[i: i + CODON_SIZE] != 'TGG':
            cleaned_sequence += seq[i: i + CODON_SIZE]
    return cleaned_sequence


def calculate_cai(prot, cu_table):
    """
    Calculate the codon adaptation index of a given protein based on a given codon usage table
    :param prot: a protein object
    :param cu_table: the codon usage table (dictionary)
    :return: the codon adaptation index
    """

    rscu_table = cu_table.get_rscu_table()
    # clean the protein from unnecessary codons
    cleaned_prot = Protein(prot.get_name(), remove_single_codons(prot.get_dna_seq()), "dna")

    # compute the cai of the actual sequence
    freq_obs = []
    for i in range(0, cleaned_prot.nucleotides_count, CODON_SIZE):
        if rscu_table[cleaned_prot.get_dna_seq()[i: i + CODON_SIZE]] > 0:
            freq_obs.append(rscu_table[cleaned_prot.get_dna_seq()[i: i + CODON_SIZE]])
        else:
            freq_obs.append(0.5)
    cai_obs = math.pow(math.prod(freq_obs), 1 / cleaned_prot.total_amino_acids_count)

    # compute the cai based on the assumption of the use of the preferred codon for each AA
    freq_max = []
    preferred_codon_rscu_for_each_aa = cu_table.get_preferred_codon_rscu_for_each_aa()
    for i in range(0, cleaned_prot.total_amino_acids_count, 1):
        freq_max.append(preferred_codon_rscu_for_each_aa[cleaned_prot.get_aa_seq()[i]])
    cai_max = math.pow(math.prod(freq_max), 1 / cleaned_prot.total_amino_acids_count)

    return round(cai_obs / cai_max,3)


if __name__ == '__main__':
    cu_table = CodonUsage("test_files/VeryHighLevelExpressionEColi.txt")
    rpsu_prot = Protein('rspu',
                        'ATGCCGGTAATTAAAGTACGTGAAAACGAGCCGTTCGACGTAGCTCTGCGTCGCTTCAAGCGTTCCTGCGAAAAAGCAGGTGTTCTGGCGGAAGTTCGTCGTCGTGAGTTCTATGAAAAACCGACTACCGAACGTAAGCGCGCTAAAGCTTCTGCAGTGAAACGTCACGCGAAGAAACTGGCTCGCGAAAACGCACGCCGCACTCGTCTGTACTAA',
                        'dna')
    print(calculate_cai(rpsu_prot, cu_table))

from utils.protein import Protein
from codon_usage.codon_usage import CodonUsage


def optimize(protein: Protein, cu_table: CodonUsage):
    """
    :param protein: the protein to optimize
    :param cu_table: the codon usage table of the species in which the protein is expressed
    :return: a fully optimize protein in which each AA is encoded by the most preferred codon
    """
    preferred_codon_for_each_aa = cu_table.get_preferred_codon_for_each_aa()

    optimized_seq = ''
    for aa in protein.get_aa_seq():
        optimized_seq += preferred_codon_for_each_aa[aa]

    optimized_prot = Protein(protein.get_name(), optimized_seq, 'DNA')

    return optimized_prot

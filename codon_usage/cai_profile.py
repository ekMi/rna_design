from codon_usage.cai import calculate_cai
from utils.protein import Protein
from codon_usage.codon_usage import CodonUsage


def cai_profile(protein: Protein, codon_window_size: int, cu_table: CodonUsage, sliding_size=1):
    """
    :param protein: the protein on which to compute the CAI profile
    :param codon_window_size: the window size on which calculate the CAI
    :param cu_table: the codon usage table of the species in which the protein is expressed
    :param sliding_size: the codon sliding size to move the window size
    :return: a list of CAI computed along the strand
    """
    cai_prof = []
    for i in range(0, len(protein.get_codon_list()) - codon_window_size, sliding_size):
        local_seq = Protein("local_seq", "".join(protein.get_codon_list()[i:i + codon_window_size]), 'dna')
        window_cai = round(calculate_cai(local_seq, cu_table), 3)
        cai_prof.append(window_cai)
    return cai_prof

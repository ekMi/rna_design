from cai import calculate_cai
from protein import Protein


def cai_profile(protein, codon_window_size, cu_table, sliding_size=1):
    cai_prof = []
    for i in range(0, len(protein.get_codon_list()) - codon_window_size, sliding_size):
        local_seq = Protein("local_seq", "".join(protein.get_codon_list()[i:i + codon_window_size]), 'dna')
        window_cai = round(calculate_cai(local_seq, cu_table), 3)
        cai_prof.append(window_cai)
    return cai_prof

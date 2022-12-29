from protein import Protein


def optimize(prot, cu_table):
    preferred_codon_for_each_aa = cu_table.get_preferred_codon_for_each_aa()

    optimized_seq = ''
    for aa in prot.get_aa_seq():
        optimized_seq += preferred_codon_for_each_aa[aa]

    optimized_prot = Protein(prot.get_name(), optimized_seq, 'DNA')

    return optimized_prot

from codon_usage import CodonUsage
from protein import Protein
from collections import namedtuple

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

CANDIDATE_AA = ['Y', 'H', 'W', 'I', 'L', 'V', 'S', 'T', 'P', 'C']


def remove_isolated_flags(flags):
    for i, flag in enumerate(flags):
        if flag:
            if i < 14:
                start = 0
            else:
                start = i - 14

            if i < len(flags) - 14:
                stop = i + 14
            else:
                stop = len(flags) - 1

            other_flag_in_region = False
            j = start
            while not other_flag_in_region and j <= stop:
                if j == i:
                    j += 1
                    continue
                else:
                    if flags[j]:
                        other_flag_in_region = True
                j += 1

            if not other_flag_in_region:
                flags[i] = False


def flag_rare_codon_for_candidate_aa(prot_org_codon_freq, threshold):
    flags = []
    # first identify codons with frequency bellow threshold AND coding for one of the CANDIDATE_AA
    for codon in prot_org_codon_freq:
        if CODON_TO_AA_DICT[codon[0]] in CANDIDATE_AA and codon[1] < threshold:
            flags.append(True)
        else:
            flags.append(False)
    print(flags)
    print(get_number_of_clusters(flags))
    remove_isolated_flags(flags)
    return flags


def get_number_of_clusters(flags):
    n_cluster = 0
    last_true_pos = -1
    for i, flag in enumerate(flags):
        if flag:
            if n_cluster == 0:
                n_cluster = 1
            # if distance from last true pos is > 14, a new cluster is starting
            elif i - last_true_pos > 14:
                n_cluster += 1
            last_true_pos = i
    return n_cluster


def a_harmonization(prot: Protein, cu_origin: CodonUsage, cu_destination: CodonUsage, threshold_frequency=0.15):
    target_number_of_clusters = round(prot.total_amino_acids_count / 30)

    cu_freq_origin = cu_origin.get_relative_codon_freq_table()

    prot_org_codon_freq = []
    for codon in prot.get_codon_list():
        prot_org_codon_freq.append((codon, cu_freq_origin[codon]))

    codon_to_modify = flag_rare_codon_for_candidate_aa(prot_org_codon_freq, threshold_frequency)



if __name__ == '__main__':
    cu_ecoli = CodonUsage("test_files/EcolCUB.txt", 'cub')
    cu_plasmodium = CodonUsage("test_files/CUB_Hive_Cut_PFaclciparum_3D7.txt", 'cub')
    protein_file = open('test_files/test_proteins.txt', 'r')

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

    msp1 = proteins[-1]
    a_harmonization(msp1, cu_plasmodium, cu_ecoli)

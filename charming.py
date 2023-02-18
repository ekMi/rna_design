from codon_rodriguez_harmonization import rodriguez_harmonization
from min_max import calculate_min_max, correlation_min_max, delta_min_max
from cai_profile import cai_profile
from codon_usage import CodonUsage
from protein import Protein
import random
from codon_randomization import randomize
from scipy.stats import pearsonr

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

CANDIDATE_SECTION_LENGTH = 5
CDR_LENGTH = 2 * CANDIDATE_SECTION_LENGTH


def calculate_distance(target_values, harmonized_values):
    total_dist = 0
    for i in range(len(target_values)):
        total_dist += abs(harmonized_values[i] - target_values[i])
    return total_dist


def calculate_above_below(target_values, harmonized_values):
    above_below = []
    for i in range(len(target_values)):
        if target_values[i] > harmonized_values[i]:
            above_below.append(0)
        elif target_values[i] < harmonized_values[i]:
            above_below.append(1)
        else:
            above_below.append(2)
    return above_below


def get_replace_pos(numReplace, replaceStart, consecutive, iteration):
    replacePos = []
    for j in range(numReplace):
        replacePos.append(
            replaceStart + int(consecutive / (numReplace + 1)) * (j + 1) - 2 + (iteration % 5))
    return replacePos


def get_replace_codon(codons_possibilities, current_codon_freq, dest_freq_table, direction):
    current_choice = None
    for i in codons_possibilities:
        if direction == 0:
            if dest_freq_table[i] > current_codon_freq:
                if current_choice is None:
                    current_choice = i
                elif dest_freq_table[i] < dest_freq_table[current_choice]:
                    current_choice = i
        else:
            if dest_freq_table[i] < current_codon_freq:
                if current_choice is None:
                    current_choice = i
                elif dest_freq_table[i] > dest_freq_table[current_choice]:
                    current_choice = i

    return current_choice

    # if direction == 0:
    #     target_rank = current_codon_rank - 1
    # else:
    #     target_rank = current_codon_rank + 1
    # for codon in codons_possibilities:
    #     if dest_ranking_table[codon] == target_rank:
    #         replace_codon = codon
    #
    # return replace_codon


def distance_to_target_reduction(protein_to_harmonize: Protein, target_vals, init_vals, dest_cu_table: CodonUsage,
                                 codon_profile_method, sliding_window):
    iteration = 0
    flag = 0
    protein_harmonized = Protein(protein_to_harmonize.get_name(), str(protein_to_harmonize.get_dna_seq()), 'dna')
    protein_harmonized_codon_list = list(protein_harmonized.get_codon_list())
    dest_freq_table = dest_cu_table.get_relative_codon_freq_table()

    while flag < CANDIDATE_SECTION_LENGTH:
        total_dist = calculate_distance(target_vals, init_vals)
        above_below = calculate_above_below(target_vals, init_vals)
        start_dist = total_dist

        # search for consecutive above_bellow_vals
        last = -1
        consecutive = 0
        for i in range(len(above_below)):
            if above_below[i] == last:
                consecutive += 1
            else:
                if consecutive >= CDR_LENGTH and last != 2:
                    prot_seq_harmo = list(protein_harmonized_codon_list)
                    num_replace = int(consecutive / CDR_LENGTH)
                    replace_start = i - consecutive
                    replace_pos = get_replace_pos(num_replace, replace_start, consecutive, iteration)
                    for j in replace_pos:
                        codon_to_replace = protein_harmonized_codon_list[j]
                        corresponding_AA = CODON_TO_AA_DICT[codon_to_replace]
                        codons_possibilities = AA_TO_CODON_DICT[corresponding_AA]
                        codon_to_replace_freq = dest_freq_table[codon_to_replace]
                        direction = above_below[i - 1]
                        replace_codon = get_replace_codon(codons_possibilities, codon_to_replace_freq,
                                                          dest_freq_table, direction)
                        if replace_codon is not None:
                            prot_seq_harmo[j] = replace_codon

                        new_prot = Protein(protein_harmonized.get_name(), "".join(prot_seq_harmo), 'dna')
                        harmonized_values = codon_profile_method(new_prot, sliding_window, dest_cu_table)
                        harmonized_dist = calculate_distance(target_vals, harmonized_values)
                        if harmonized_dist < total_dist:
                            protein_harmonized = new_prot
                            protein_harmonized_codon_list = new_prot.get_codon_list()
                            init_vals = list(harmonized_values)
                            total_dist = harmonized_dist
                            flag = 0
                last = above_below[i]
                consecutive = 1
        if start_dist == total_dist:
            flag += 1
        iteration += 1
    return protein_harmonized


def generate_random_seq(input_protein: Protein, nb_of_sequences):
    random_protein_sequences = []
    input_aa_seq = input_protein.get_aa_seq()
    for i in range(nb_of_sequences):
        seq = []
        for aa in input_aa_seq:
            seq.append(AA_TO_CODON_DICT[aa][random.randint(0, len(AA_TO_CODON_DICT[aa]) - 1)])
        prot = Protein(f'{input_protein.get_name()}Random{str(i)}', "".join(seq), 'dna')
        random_protein_sequences.append(prot)

    return random_protein_sequences


def generate_random_seq_freq_conserved(input_protein: Protein, nb_of_sequences: int, dest_cu: CodonUsage):
    random_protein_sequences = []
    for i in range(nb_of_sequences):
        prot = randomize(input_protein, dest_cu)
        prot.set_name(f'{prot.get_name()}Random{str(i)}')
        random_protein_sequences.append(prot)

    return random_protein_sequences


def charming(input_protein: Protein, init_cu_table: CodonUsage, dest_cu_table: CodonUsage, codon_profile_method='MM',
             number_output=3, sliding_window=10):
    profile_methods = {'MM': calculate_min_max, 'CAI': cai_profile}

    target_profile = profile_methods[codon_profile_method](input_protein, sliding_window, init_cu_table)

    random_protein_sequences = [Protein(f'{input_protein.get_name()}RodriguezInit', input_protein.get_dna_seq(), 'dna')]
    #
    # random_protein_sequences.extend(generate_random_seq(input_protein, number_output * 10 - 1))

    random_protein_sequences.extend(
        generate_random_seq_freq_conserved(input_protein, number_output * 10 - 1, dest_cu_table))

    harmonized_proteins = []
    for prot_id, random_prot in enumerate(random_protein_sequences):
        # use input sequence for rodriguez initialization
        if prot_id == 0:
            init_prot = rodriguez_harmonization(random_prot, init_cu_table, dest_cu_table)
            init_vals = profile_methods[codon_profile_method](init_prot, sliding_window, dest_cu_table)
        else:
            init_prot = random_prot
            init_vals = profile_methods[codon_profile_method](init_prot, sliding_window, dest_cu_table)

        # this is where the magic happens
        harmonized_prot = distance_to_target_reduction(init_prot, target_profile, init_vals, dest_cu_table,
                                                       profile_methods[codon_profile_method], sliding_window)

        final_vals = profile_methods[codon_profile_method](harmonized_prot, sliding_window, dest_cu_table)
        corr, _ = pearsonr(target_profile, final_vals)
        harmonized_proteins.append((corr, harmonized_prot))

    harmonized_proteins.sort(reverse=True)

    output_proteins = []
    for i in range(number_output):
        output_proteins.append(harmonized_proteins[i][1])

    return output_proteins

from codon_usage import CodonUsage
from protein import Protein
from cai import calculate_cai
from min_max import calculate_min_max
import pandas as pd
import plotly.express as px

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

def rodriguez_harmonization(protein, origin_cu_table, destination_cu_table):

    """
    Harmonize the sequence to the destiantion host based on the codon ranking for each amino acid
    :param protein: the protein to harmonize
    :param origin_cu_table: the cu table of the inital host
    :param destination_cu_table: the cu table of the detination host
    :return: an harmonized protein to the destination host based on the codon ranking of the inital host
    """

    print(protein.get_name())

    codon_list = protein.get_codon_list()
    origin_codon_ranking_table = origin_cu_table.get_relative_codon_ranking_table()
    dest_codon_ranking_table = destination_cu_table.get_relative_codon_ranking_table()

    harmonized_sequence = ""
    # Reads the initial sequence codon by codon
    for i, init_codon in enumerate(codon_list):
        # Get the other codons coding for the same AA
        aa = CODON_TO_AA_DICT[init_codon]
        possible_codons = AA_TO_CODON_DICT[aa]
        # Loop over those possibilities
        for codon in possible_codons:
            # Find the one for which the rank considering the destination cu table
            # is equal to the rank of the initial codon with the origin cu table
            if dest_codon_ranking_table[codon] == origin_codon_ranking_table[init_codon]:
                # Append it to the harmonized sequence
                harmonized_sequence += codon

    # Use the harmonized sequence to create a protein object
    harmonized_prot = Protein(protein.get_name(), harmonized_sequence, 'dna')

    return harmonized_prot

if __name__ == '__main__':
    cu_ecoli = CodonUsage("test_files/EcolCUB.txt", 'cub')
    cu_human = CodonUsage("test_files/CUB_HiveCut_HSapiens.txt", 'cub')

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


    alpha_actin = proteins[- 2]
    alpha_actin = Protein(alpha_actin.get_name(), "".join(alpha_actin.get_codon_list()[0:250]), 'dna')
    alpha_actin_harmonized = rodriguez_harmonization(alpha_actin, cu_human, cu_ecoli)

    min_max_wt_in_human = calculate_min_max(alpha_actin, 17, cu_human)
    min_max_wt_in_ecoli = calculate_min_max(alpha_actin, 17, cu_ecoli)
    min_max_harmonized_ecoli = calculate_min_max(alpha_actin_harmonized, 17, cu_ecoli)

    human_rank = cu_human.get_relative_codon_ranking_table()
    ecoli_rank = cu_ecoli.get_relative_codon_ranking_table()


    df = pd.DataFrame(zip(list(min_max_wt_in_human), list(min_max_wt_in_ecoli), list(min_max_harmonized_ecoli)) , columns=['wt_human', 'wt_ecoli', 'harmo_ecoli'])

    fig = px.line(df)
    fig.update_yaxes(range=[-100, 100])
    fig.show()

    print(calculate_cai(alpha_actin, cu_human))
    print(calculate_cai(alpha_actin_harmonized, cu_human))
    print(calculate_cai(alpha_actin_harmonized, cu_ecoli))

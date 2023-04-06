from utils.codon_aa_dict import AA_TO_CODON_DICT, CODON_TO_AA_DICT
from codon_usage.codon_usage import CodonUsage
from utils.protein import Protein


def rodriguez_harmonization(protein: Protein, origin_cu_table: CodonUsage, destination_cu_table: CodonUsage):
    """
    Harmonize the sequence to the destination host based on the codon ranking for each amino acid
    :param protein: the protein to harmonize
    :param origin_cu_table: the cu table of the initial host
    :param destination_cu_table: the cu table of the destination host
    :return: an harmonized protein to the destination host based on the codon ranking of the initial host
    """

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

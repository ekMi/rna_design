CODON_SIZE = 3

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
AA_LIST = ['S', 'N', 'W', 'D', 'E', 'F', 'P', 'T', 'G', 'R', 'L', 'V', 'A', '*', 'H', 'I', 'K', 'Y', 'M', 'C', 'Q']


def convert_dna_to_rna(sequence):
    return sequence.replace('T', 'U')


def convert_rna_to_dna(sequence):
    return sequence.replace('U', 'T')


def convert_dna_to_amino_acids(dna_seq):
    amino_acid_sequence = ''
    for i in range(0, len(dna_seq), CODON_SIZE):
        amino_acid_sequence += CODON_TO_AA_DICT[dna_seq[i:i + CODON_SIZE]]
    return amino_acid_sequence


class Protein:
    """
    Class containing protein information:
    - DNA coding sequence
    - RNA coding sequence
    - AA sequence
    - Number of nucleotides
    - Number of total amino acids
    - Number of each amino acid (among 20 possibilities)

    Can be instantiated with a DNA or RNA sequence and its name
    """

    def __init__(self, name: str, sequence: str, sequence_type: str):
        if sequence_type.upper() == 'DNA':
            if len(sequence) % CODON_SIZE == 0 and all(x in ['A', 'T', 'G', 'C'] for x in sequence):
                self._DNA_sequence = sequence.upper()
                self._RNA_sequence = convert_dna_to_rna(self._DNA_sequence)
                self._AA_sequence = convert_dna_to_amino_acids(self._DNA_sequence)
            else:
                raise ValueError("Protein DNA sequence must be of a multiple of 3 and contain only A-T-G-C letters")
        elif sequence_type.upper() == 'RNA' and all(x in ['A', 'U', 'G', 'C'] for x in sequence):
            if len(sequence) % CODON_SIZE == 0:
                self._RNA_sequence = sequence
                self._DNA_sequence = convert_rna_to_dna(self._RNA_sequence)
                self._AA_sequence = convert_dna_to_amino_acids(self._DNA_sequence)
            else:
                raise ValueError("Protein nucleotide sequence must be of a multiple of 3")
        else:
            raise ValueError('Protein sequence_type must be DNA, RNA')
        self._name = name

    @property
    def nucleotides_count(self):
        return len(self._DNA_sequence)

    @property
    def total_amino_acids_count(self):
        return len(self._AA_sequence)

    @property
    def amino_acid_counts(self):
        '''
        :return: A dictionary with the number of occurence for each amino acid
        '''
        aa_count_dict = {'S': 0, 'N': 0, 'W': 0, 'D': 0, 'E': 0, 'F': 0, 'P': 0, 'T': 0, 'G': 0, 'R': 0, 'L': 0,
                         'V': 0, 'A': 0, '*': 0, 'H': 0, 'I': 0, 'K': 0, 'Y': 0, 'M': 0, 'C': 0, 'Q': 0}

        for aa, value in aa_count_dict.items():
            aa_count_dict[aa] = self._AA_sequence.count(aa)

        return aa_count_dict

    def get_codon_list(self):
        codon_list = []
        for i in range(0, len(self._DNA_sequence), CODON_SIZE):
            codon_list.append(self._DNA_sequence[i:i+CODON_SIZE])
        return codon_list

    def get_rna_seq(self):
        return self._RNA_sequence

    def get_dna_seq(self):
        return self._DNA_sequence

    def get_aa_seq(self):
        return self._AA_sequence

    def get_name(self):
        return self._name



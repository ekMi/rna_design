from utils.codon_aa_dict import AA_TO_CODON_DICT, CODON_TO_AA_DICT, CODON_SIZE

CODON_LIST = ['TCA', 'AAT', 'TGG', 'GAT', 'GAA', 'TTC', 'CCG', 'ACT', 'GGG', 'ACG', 'AGA', 'TTG', 'GTC', 'GCA',
              'TGA', 'CGT', 'CAC', 'CTC', 'CGA', 'GCT', 'ATC', 'ATA', 'TTT', 'TAA', 'GTG', 'GCC', 'GAG', 'CAT',
              'AAG', 'AAA', 'GCG', 'TCC', 'GGC', 'TCT', 'CCT', 'GTA', 'AGG', 'CCA', 'TAT', 'ACC', 'TCG', 'ATG',
              'TTA', 'TGC', 'GTT', 'CTT', 'CAG', 'CCC', 'ATT', 'ACA', 'AAC', 'GGT', 'AGC', 'CGG', 'TAG', 'CGC',
              'AGT', 'CTA', 'CAA', 'CTG', 'GGA', 'TGT', 'TAC', 'GAC']


class CodonUsage:
    """
    Compute codon usage table based on a reference sequence contained in a *.txt file (with DNA nt only)
    or from a cub table retrieved from a database
    This table serve as the basis to compute the RSCU, relative adaptiveness, relative codon frequency tables,
    relative codon ranking table
    """

    def __init__(self, file, file_format="sequence"):

        full_file = open(file, 'r')
        lines = full_file.readlines()

        if file_format == "sequence":
            # The cub table must be computed from the full string of the file
            genome = ""
            for line in lines:
                if not ":" in line:  # skip line containing protein name
                    genome += line.strip()

            nucleotide_count = len(genome)
            if nucleotide_count % CODON_SIZE != 0 or not all(x in ['A', 'T', 'G', 'C'] for x in genome):
                raise ValueError('calculate_codon_usage_table: Sequence must be of a length multiple of 3 and contain '
                                 'only ATGC letters')

            aa_count = nucleotide_count / CODON_SIZE

            self._codon_usage_table = dict((k, 0) for k in CODON_LIST)
            for i in range(0, nucleotide_count, CODON_SIZE):
                self._codon_usage_table[genome[i: i + CODON_SIZE]] += 1

            self._codon_usage_table = dict(
                (k, round((v / aa_count) * 1000, 2)) for k, v in self._codon_usage_table.items())

        elif file_format == "cub":
            self._codon_usage_table = dict()
            for line in lines:
                line = line.split()
                if len(line) < 2:
                    raise ValueError("calculate_codon_usage_table: Invalid record format, must be codon name value "
                                     "pair separated by a space")
                if line[0] in CODON_TO_AA_DICT:
                    self._codon_usage_table[line[0]] = float(line[1])

            if len(self._codon_usage_table) < len(CODON_TO_AA_DICT):
                missing_count = len(CODON_TO_AA_DICT) - len(self._codon_usage_table)
                raise ValueError(f"calculate_codon_usage_table: There are {missing_count} missing codons in the file "
                                 f"provided")
        else:
            raise ValueError("calculate_codon_usage_table: The 2nd parameter for file format must be 'sequence' or "
                             "'cub' depending on the input file")

    def get_codon_usage_table(self):
        """
        :return: the computed codon usage table
        """
        return self._codon_usage_table

    def get_rscu_table(self):
        """
        Compute the rscu table based on the codon usage table. For each codon, RSCU = codon usage value / (sum(values
        for each codon coding for the same aa)/number of codon coding for the same aa) :return: RSCU_table
        """
        aa_sum_freq = {'S': 0, 'N': 0, 'W': 0, 'D': 0, 'E': 0, 'F': 0, 'P': 0, 'T': 0, 'G': 0, 'R': 0, 'L': 0,
                       'V': 0, 'A': 0, '*': 0, 'H': 0, 'I': 0, 'K': 0, 'Y': 0, 'M': 0, 'C': 0, 'Q': 0}

        for key, value in aa_sum_freq.items():
            aa_sum_freq[key] = sum([v for k, v in self._codon_usage_table.items() if k in AA_TO_CODON_DICT[key]])

        return dict((codon_key, round(codon_value / (sum_codon_freq_for_aa / len(codons)), 3))
                    for codon_key, codon_value in self._codon_usage_table.items()
                    for aa, sum_codon_freq_for_aa in aa_sum_freq.items() if aa == CODON_TO_AA_DICT[codon_key]
                    for aa2, codons in AA_TO_CODON_DICT.items() if aa2 == CODON_TO_AA_DICT[codon_key])

    def get_relative_adaptiveness_table(self):
        """
        Compute the relative adaptiveness table based on the codon usage table. The relative adaptiveness  for each
        codon = RSCU codon / RSCU max of codon coding for same aa = codon usage value / max codon usage value for the
        same aa :return: relative adaptiveness table as a disctionary
        """
        return dict((codon_key, round(codon_value / max_codon, 3))
                    for codon_key, codon_value in self._codon_usage_table.items()
                    for max_codon in [max([v for key, v in self._codon_usage_table.items() if
                                           key in AA_TO_CODON_DICT[CODON_TO_AA_DICT[codon_key]]])])

    def get_preferred_codon_rscu_for_each_aa(self):
        """
        :return: a dictionary with the preferred codon for each aa based on the rscu values
        """
        rscu_table = self.get_rscu_table()
        preferred_codon_for_each_aa = {}
        for aa, _ in AA_TO_CODON_DICT.items():
            preferred_codon_for_each_aa[aa] = max([v for k, v in rscu_table.items() if k in AA_TO_CODON_DICT[aa]])
        return preferred_codon_for_each_aa

    def get_preferred_codon_for_each_aa(self):
        """
        :return: a dictionary with the preferred codon for each aa based on the cu values
        """
        preferred_codon_for_each_aa = {}
        for aa, _ in AA_TO_CODON_DICT.items():
            values = dict((k, v) for k, v in self._codon_usage_table.items() if k in AA_TO_CODON_DICT[aa])
            preferred_codon_for_each_aa[aa] = max(values, key=values.get)
        return preferred_codon_for_each_aa

    def get_relative_codon_freq_table(self):
        """
        Computes the relative codon frequency table based on the rscu table. The relative codon frequency is the
        computed as the RSCU value of a codon divided by the number of possibilities of codon for its corresponding AA
        :return: a dictionnary with the relative frequency of each codon
        """
        rscu_table = self.get_rscu_table()
        codon_freq_dict = dict(
            (key, round(value / len(codons_for_aa), 2)) for key, value in rscu_table.items() for aa, codons_for_aa in
            AA_TO_CODON_DICT.items() if key in codons_for_aa)
        return codon_freq_dict

    def get_relative_codon_ranking_table(self):
        """
        Calculate the codon preference (from 1 to 6) order relative to its corresponding AA
        :return: a dictionary with the ranking for each codon
        """

        codon_ranking_dict = dict((key, 0) for key in CODON_TO_AA_DICT)

        for aa, codons in AA_TO_CODON_DICT.items():
            codon_val_dict = dict((codon, self._codon_usage_table[codon]) for codon in codons)
            rank = 1
            for codon, value in sorted(codon_val_dict.items(), key=lambda v: v[1], reverse=True):
                codon_ranking_dict[codon] = rank
                rank += 1
        return codon_ranking_dict


if __name__ == '__main__':
    cu_table = CodonUsage("../test/test_files/VeryHighLevelExpressionEColi.txt")
    print(sorted(cu_table.get_rscu_table().items()))

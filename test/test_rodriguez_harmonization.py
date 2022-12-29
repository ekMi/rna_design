from protein import Protein
from codon_usage import CodonUsage
from cai import calculate_cai
from min_max import calculate_min_max, delta_min_max, correlation_min_max
from codon_rodriguez_harmonization import rodriguez_harmonization
import plotly.express as px
import pandas as pd

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

if __name__ == '__main__':
    CODON_WINDOW_SIZE = 10
    CODON_SHIFT = 1

    cu_table_ecoli = CodonUsage('../test_files/CUB_HiveCut_Ecoli_K12.txt', 'cub')
    cu_table_scer = CodonUsage('../test_files/ScerCUB.txt', 'cub')
    protein_file = open('../test_files/test_proteins.txt', 'r')

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

    # test on the third protein (rpoD) which is relatively long, so use only the 200 first AAs for graph visibility
    rpoD = proteins[2]
    rpoD = Protein(rpoD.get_name(), "".join(rpoD.get_codon_list()[0:200]), 'dna')

    # get the frequency and codon of the codon at Rank1 for Asparagine (R) in Ecoli and Scer cu table
    codons_for_arg = AA_TO_CODON_DICT['R']
    ecoli_first_codon = cu_table_ecoli.get_preferred_codon_for_each_aa()['R']
    scer_first_codon = cu_table_scer.get_preferred_codon_for_each_aa()['R']

    print(f'The codon at Rank1 for Arginine in Ecoli is: {ecoli_first_codon}, it\'s frequency is: {cu_table_ecoli.get_relative_codon_freq_table()[ecoli_first_codon]}')
    print(f'The codon at Rank1 for Arginine in Scer is: {scer_first_codon}, it\'s frequency is: {cu_table_scer.get_relative_codon_freq_table()[scer_first_codon]}')

    minmax_native_ecoli = calculate_min_max(rpoD, CODON_WINDOW_SIZE, cu_table_ecoli,1)
    minmax_native_scer = calculate_min_max(rpoD, CODON_WINDOW_SIZE, cu_table_scer, 1)

    rpoD_harmonized = rodriguez_harmonization(rpoD, cu_table_ecoli, cu_table_scer)
    min_max_harmonized_scer = calculate_min_max(rpoD_harmonized, CODON_WINDOW_SIZE, cu_table_scer, 1)
    print(delta_min_max(minmax_native_ecoli, minmax_native_scer))
    print(delta_min_max(minmax_native_ecoli, min_max_harmonized_scer))
    print(correlation_min_max(minmax_native_ecoli, minmax_native_scer))
    print(correlation_min_max(minmax_native_ecoli, min_max_harmonized_scer))



    df = pd.DataFrame(zip(list(minmax_native_ecoli), list(minmax_native_scer), list(min_max_harmonized_scer)) , columns=['wt_ecoli', 'wt_scer', 'harmo_scer'])

    fig = px.line(df)
    fig.update_yaxes(range=[-100, 100])

    fig.write_image("images/fig_rpoD_rodriguez_scer_minmax.jpeg")


    print(calculate_cai(rpoD, cu_table_ecoli))
    print(calculate_cai(rpoD, cu_table_scer))
    print(calculate_cai(rpoD_harmonized, cu_table_scer))
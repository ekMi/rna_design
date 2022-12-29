from protein import Protein
from codon_usage import CodonUsage
from cai import calculate_cai
from min_max import calculate_min_max, delta_min_max, correlation_min_max
from charming import charming
import plotly.express as px
import pandas as pd


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

    harmonized_rpoD = charming(rpoD, cu_table_ecoli, cu_table_scer, 'MM', 3, 10)

    minmax_native_ecoli = calculate_min_max(rpoD, CODON_WINDOW_SIZE, cu_table_ecoli,1)
    minmax_native_scer = calculate_min_max(rpoD, CODON_WINDOW_SIZE, cu_table_scer, 1)

    print('native sequence/native host')
    print(calculate_cai(rpoD, cu_table_ecoli))
    print('native seequence/ destination host')
    print(calculate_cai(rpoD, cu_table_scer))
    print(delta_min_max(minmax_native_ecoli, minmax_native_scer))
    print(correlation_min_max(minmax_native_ecoli, minmax_native_scer))


    min_max_harmonized = []
    for i in range(3):
        min_max_harmo = calculate_min_max(harmonized_rpoD[i], CODON_WINDOW_SIZE, cu_table_scer, 1)
        min_max_harmonized.append(min_max_harmo)
        print(harmonized_rpoD[i].get_name())
        print(calculate_cai(harmonized_rpoD[i], cu_table_scer))
        print(delta_min_max(minmax_native_ecoli, min_max_harmo))
        print(correlation_min_max(minmax_native_ecoli, min_max_harmo))

    df = pd.DataFrame(zip(list(minmax_native_ecoli), list(minmax_native_scer), list(min_max_harmonized[0]), list(min_max_harmonized[1]), list(min_max_harmonized[2])),
                      columns=['wt_ecoli', 'wt_scer', harmonized_rpoD[0].get_name(), harmonized_rpoD[1].get_name(), harmonized_rpoD[2].get_name()])

    fig = px.line(df)
    fig.update_yaxes(range=[-100, 100])
    fig.write_image("images/fig_rpoD_charming_scer_minmax.jpeg")


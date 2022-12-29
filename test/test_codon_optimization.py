from protein import Protein
from codon_usage import CodonUsage
from cai import calculate_cai
from min_max import calculate_min_max, correlation_min_max, delta_min_max
from codon_optimization import optimize
import plotly.express as px
import pandas as pd

if __name__ == '__main__':
    cu_table = CodonUsage('../test_files/CUB_HiveCut_Ecoli_K12.txt', 'cub')
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
    proteins.append(Protein(protein_name,sequence,'dna'))

    for protein in proteins:
        prot = optimize(protein, cu_table)
        print(prot.get_name())
        print(round(calculate_cai(prot, cu_table), 3))

    rpoD = proteins[2]
    rpoD = Protein(rpoD.get_name(), "".join(rpoD.get_codon_list()[0:200]), 'dna')
    rpoD_optimized_scer = optimize(rpoD, cu_table_scer)

    minmax_native_ecoli = calculate_min_max(rpoD, 10, cu_table)
    minmax_native_scer = calculate_min_max(rpoD, 10, cu_table_scer)
    minmax_optimized_scer = calculate_min_max(rpoD_optimized_scer, 10, cu_table_scer)

    print(delta_min_max(minmax_native_ecoli, minmax_native_ecoli))
    print(delta_min_max(minmax_native_ecoli, minmax_native_scer))
    print(delta_min_max(minmax_native_ecoli, minmax_optimized_scer))
    print(correlation_min_max(minmax_native_ecoli, minmax_native_ecoli))
    print(correlation_min_max(minmax_native_ecoli, minmax_native_scer))
    print(correlation_min_max(minmax_native_ecoli, minmax_optimized_scer))

    df = pd.DataFrame(zip(list(minmax_native_ecoli), list(minmax_native_scer), list(minmax_optimized_scer)) , columns=['wt_ecoli', 'wt_scer', 'optimized_scer'])

    fig = px.line(df)
    fig.update_yaxes(range=[-100, 100])

    fig.write_image("images/fig_rpoD_optimized_scer_minmax.jpeg")


    print(calculate_cai(rpoD, cu_table))
    print(calculate_cai(rpoD, cu_table_scer))
    print(calculate_cai(rpoD_optimized_scer, cu_table_scer))



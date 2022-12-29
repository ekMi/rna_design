from protein import Protein
from codon_usage import CodonUsage
from cai import calculate_cai
from min_max import calculate_min_max, correlation_min_max, delta_min_max
from codon_randomization import randomize
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

    print(calculate_cai(rpoD, cu_table_ecoli))
    print(calculate_cai(rpoD, cu_table_scer))
    minmax_native_ecoli = calculate_min_max(rpoD, CODON_WINDOW_SIZE, cu_table_ecoli,1)
    minmax_native_scer = calculate_min_max(rpoD, CODON_WINDOW_SIZE, cu_table_scer, 1)
    print(delta_min_max(minmax_native_ecoli, minmax_native_ecoli))
    print(correlation_min_max(minmax_native_ecoli, minmax_native_ecoli))
    print(delta_min_max(minmax_native_ecoli, minmax_native_scer))
    print(correlation_min_max(minmax_native_ecoli, minmax_native_scer))




    # create 2 randomized sequence
    min_max_profiles_randomized = []
    for i in range(2):
        cai_profile_randomized = []
        randomized_rpoD = randomize(rpoD, cu_table_scer)
        print(calculate_cai(randomized_rpoD, cu_table_scer))
        random_minmax = calculate_min_max(randomized_rpoD, CODON_WINDOW_SIZE, cu_table_scer,1)
        min_max_profiles_randomized.append(random_minmax)
        print(delta_min_max(minmax_native_ecoli, random_minmax))
        print(correlation_min_max(minmax_native_ecoli, random_minmax))


    df = pd.DataFrame(list(zip(minmax_native_ecoli, minmax_native_scer ,min_max_profiles_randomized[0], min_max_profiles_randomized[1])),
                      columns=['Native sequence from E coli with E coli usage table', 'Native sequence from E coli with S cer usage table','Randomized sequence 1 with S cer usage table',
                               'Randomized sequence 2 with S cer usage table'])

    fig = px.line(df)
    fig.update_yaxes(range=[-100, 100])
    fig.update_layout(
        xaxis_title="Codon windows (10 codons /window)", yaxis_title="%MinMax"
    )
    fig['data'][0]['line']['color'] = "#000000"
    fig['data'][0]['line']['dash'] = "dash"
    fig['data'][1]['line']['color'] = "#000000"
    fig['data'][2]['line']['color'] = "#000FFF"
    fig['data'][3]['line']['color'] = "#1FFF00"
    fig.update_layout(legend=dict(
        yanchor="bottom",
        y=0.05,
        xanchor="left",
        x=0.01,
        font = dict(size = 10)
    ))


    fig.show()

    fig.write_image("images/fig_rpoD_randomized_scer_minmax.jpeg")

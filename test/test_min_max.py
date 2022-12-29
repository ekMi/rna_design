from protein import Protein
from codon_usage import CodonUsage
from min_max import calculate_min_max
from cai_profile import cai_profile
import pandas as pd
import plotly.express as px


if __name__ == '__main__':
    CODON_WINDOW_SIZE = 10
    SLIDING_SIZE = 1

    cu_table = CodonUsage('../test_files/CUB_HiveCut_Ecoli_K12.txt', 'cub')
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

    # test on the third protein (rpoD) which is relatively long
    rpoD = proteins[2]

    rpoD_min_max = calculate_min_max(rpoD, CODON_WINDOW_SIZE, cu_table, SLIDING_SIZE)

    df = pd.DataFrame( list(rpoD_min_max) , columns=['%MinMax'])

    fig = px.line(df)
    fig.update_yaxes(range=[-100, 100])
    fig.update_layout(
    xaxis_title="Codon windows (10 codons/window sliding over 1 codons)", yaxis_title="%MinMax"
    )

    fig.write_image("images/fig_rpoD_MinMax_profile_W10_S1.jpeg")

    rpoD_cai_profile = cai_profile(rpoD, cu_table, CODON_WINDOW_SIZE, SLIDING_SIZE)

    df = pd.DataFrame( list(rpoD_cai_profile) , columns=['CAI profile'])

    fig = px.line(df)
    fig.update_yaxes(range=[0, 1])
    fig.update_layout(
    xaxis_title="Codon windows (1O codons/window sliding over 1 codon)", yaxis_title="CAI profile"
    )

    fig.write_image("images/fig_rpoD_CAI_profile_W10_S1.jpeg")

    # freqDict = readCUBTable('../test_files/CUB_HiveCut_Ecoli_K12.txt')
    # aaFreqDict = generateAAFreqDict(freqDict, aaDict)
    #
    # print(f'Testing minmax calculation 100 time over {len(proteins)}')
    # print('Optimized algorithm')
    # start = time.time()
    # for protein in proteins:
    #     for i in range(1000):
    #         calculate_min_max(protein,CODON_WINDOW_SIZE, cu_table, SLIDING_SIZE)
    # print(f'Took {time.time() - start}')
    #
    #
    # print('Authors algorithm')
    # start = time.time()
    # for protein in proteins:
    #     for i in range(1000):
    #         calculateMinMax(protein.get_codon_list(),aaFreqDict, freqDict, mapDict, 18)
    # print(f'Took {time.time() - start}')



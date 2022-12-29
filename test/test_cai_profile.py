from protein import Protein
from codon_usage import CodonUsage
from cai import calculate_cai
import plotly.express as px
import pandas as pd

if __name__ == '__main__':
    CODON_WINDOW_SIZE = 15
    SLIDING_SIZE = 5

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
    cai_profile = []
    cai_global = []
    global_cai = round(calculate_cai(rpoD, cu_table), 3)
    print(global_cai)
    print(rpoD.get_name())
    for i in range(0, len(rpoD.get_codon_list()) - CODON_WINDOW_SIZE, SLIDING_SIZE):
        local_seq = Protein("local_seq", "".join(rpoD.get_codon_list()[i:i + CODON_WINDOW_SIZE]), 'dna')
        window_cai = round(calculate_cai(local_seq, cu_table), 3)
        cai_profile.append(window_cai)
        cai_global.append(global_cai)

    df = pd.DataFrame( list(zip(cai_profile, cai_global) ), columns=['CAI Profile', 'Global CAI'])

    fig = px.line(df)
    fig.update_yaxes(range=[0, 1])
    fig.update_layout(
    xaxis_title="Codon windows (15 codons /window)", yaxis_title="CAI"
    )
    # fig.show()

    # fig.write_image("images/fig_rpoD_CAI_profile.jpeg")

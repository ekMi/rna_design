from codon_usage import CodonUsage
import plotly.express as px
import pandas as pd
import os

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
    cub_hivecut = CodonUsage('../test_files/CUB_HiveCut_Ecoli_K12.txt', 'cub').get_relative_codon_freq_table()
    cub_kazusa = CodonUsage('../test_files/CUB_Kazusa_Ecoli_K12.txt', 'cub').get_relative_codon_freq_table()
    cub_computed = CodonUsage('../test_files/VeryHighLevelExpressionEColi.txt').get_relative_codon_freq_table()

    heatmap = pd.DataFrame(columns=['Amino Acid', 'Codon', 'Kazusa', 'Hive Cut', 'Computed'])

    for aa in AA_TO_CODON_DICT:
        for codon in AA_TO_CODON_DICT[aa]:
            df = pd.DataFrame([[aa, codon, cub_kazusa[codon], cub_hivecut[codon], cub_computed[codon]]],
                              columns=['Amino Acid', 'Codon', 'Kazusa', 'Hive Cut', 'Computed'])
            heatmap = pd.concat([heatmap, df])
    heatmap.reset_index(inplace=True, drop=True)
    data = heatmap[['Kazusa', 'Hive Cut', 'Computed']]

    fig1 = px.imshow(data.head(34), y=heatmap['Amino Acid'].head(34) + " " + heatmap['Codon'].head(34), color_continuous_scale="Greens")
    fig1.layout.yaxis.type = "category"
    fig1.layout.autosize = True
    fig1.layout.yaxis.tickmode = 'linear'
    fig1.layout.xaxis.tickmode = 'linear'
    fig1.layout.xaxis.tickangle = 65
    fig1.update_layout(
        yaxis=dict(
            tickfont=dict(size=10)))
    fig1.update_layout(
        xaxis=dict(
            tickfont=dict(size=10)))



    fig2 = px.imshow(data.tail(30), y=heatmap['Amino Acid'].tail(30) + " " + heatmap['Codon'].tail(30), color_continuous_scale="Greens")
    fig2.layout.yaxis.type = "category"
    fig2.layout.autosize = True
    fig2.layout.yaxis.tickmode = 'linear'
    fig2.layout.xaxis.tickmode = 'linear'
    fig2.layout.xaxis.tickangle = 65
    fig2.update_layout(
        yaxis=dict(
            tickfont=dict(size=10)))
    fig2.update_layout(
        xaxis=dict(
            tickfont=dict(size=10)))

    if not os.path.exists("images"):
        os.mkdir("images")

    fig1.write_image("images/fig1.jpeg")
    fig2.write_image("images/fig2.jpeg")


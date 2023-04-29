from utils.protein import Protein
from utils.import_proteins_file import import_proteins
from codon_usage.codon_usage import CodonUsage
from codon_ramp import replace_codons_with_rare
from codon_usage.min_max import calculate_min_max, delta_min_max, correlation_min_max
from codon_usage.cai import calculate_cai
from codon_bias_optimizations.charming import charming
from structure_reduction.struct_optimizations import SimulatedAnnealing, StructReduction
import pandas as pd
import plotly.express as px
import RNA
from alternative_start_removal import find_alternative_starts, replace_codons
from PIL import Image


if __name__ == '__main__':
    CODON_WINDOW_SIZE = 10
    CODON_SHIFT = 1

    # loading codon usage tables
    cu_table_human = CodonUsage('test_files/CUB_HiveCut_HSapiens.txt', 'cub')       # the origin CU table
    cu_table_ecoli = CodonUsage('test_files/CUB_HiveCut_Ecoli_K12.txt', 'cub')                      # the destination CU table

    # get proteins (encoded with DNA string)
    proteins = import_proteins("test_files/hepatitis_e_virus.fasta")
    orf2 = proteins[0]

    # split the seq in two parts
    # the first 48 nucleotides to lower the codon usage frequency
    orf2_start_low_codons_freq = replace_codons_with_rare(orf2.get_dna_seq()[:48], cu_table_ecoli)

    print("Sequence 16 first codons freq lowered")

    # the rest to harmonize with dest cu_table
    orf2_rest = Protein(orf2.get_name(), ''.join(orf2.get_codon_list()[16:]), 'dna')
    orf2_rest_harmonized = charming(orf2_rest, cu_table_human, cu_table_ecoli, 'MM', 1, CODON_WINDOW_SIZE)

    print("Sequence harmonized")

    # recombine in one single protein
    orf2_full_harmonized = Protein(
        orf2.get_name(),
        orf2_start_low_codons_freq + orf2_rest_harmonized[0].get_dna_seq(),
        'dna'
    )

    # remove alternative starts
    alt_starts = find_alternative_starts(orf2_full_harmonized.get_dna_seq())
    new_seq = replace_codons(orf2_full_harmonized.get_dna_seq(), alt_starts)

    print("Alternative start removed")

    orf2_alt_start_removed = Protein(
        orf2.get_name(),
        new_seq,
        'dna'
    )

    struct_reduction_strategy = SimulatedAnnealing(max_iteration=1000)
    struct_reduction = StructReduction(struct_reduction_strategy, cu_table_ecoli)

    sd = 'AGGAGG'
    spacer = 'AAAAAAAA'

    best_spacer, best_seq, best_score, iter_best = struct_reduction.reduce_structure(orf2_alt_start_removed.get_dna_seq(), sd, spacer, 0.4)
    print(f'Best score {best_score}')
    print(f'found after {iter_best} iterations')
    structure, mfe = RNA.fold(sd + best_spacer + best_seq)
    print('Best sequence:')
    print(sd + best_spacer + best_seq)
    print('Best structure:')
    print(structure)


    # Génération de la notation PostScript
    ps_file = "orf2.ps"
    RNA.PS_rna_plot(sd + best_spacer + best_seq, structure, ps_file)

    # Convertir la notation PostScript en image PNG
    img_file = "rna_folding_optimized_orf2.png"

    img = Image.open(ps_file)
    img.save(img_file, "png")

    best_seq = Protein(orf2.get_name(), best_seq, 'dna')

    # calculate minmax
    minmax_native_human = calculate_min_max(orf2, CODON_WINDOW_SIZE, cu_table_human, CODON_SHIFT)
    minmax_native_ecoli = calculate_min_max(orf2, CODON_WINDOW_SIZE, cu_table_ecoli, CODON_SHIFT)
    min_max_full_best_seq = calculate_min_max(best_seq, CODON_WINDOW_SIZE, cu_table_ecoli, CODON_SHIFT)

    # output the result in a chart
    df = pd.DataFrame(zip(list(minmax_native_human), list(minmax_native_ecoli), list(min_max_full_best_seq),),
                      columns=['wt_human', 'wt_ecoli', 'optimized_ecoli'])
    # compute the ratio
    df['reference'] = df['wt_human'] - df['wt_human']
    df['ratio_wild_type_ecoli'] = df['wt_ecoli'] - df['wt_human']
    df['ratio_optimized_ecoli'] = df['optimized_ecoli'] - df['wt_human']

    print(df)

    df.drop(['wt_human', 'wt_ecoli', 'optimized_ecoli'], inplace=True, axis=1)

    fig = px.line(df)
    fig.update_yaxes(range=[-100, 100])
    fig.write_image("fig_orf2_full_workflow.jpeg")
    fig.show()

    # compute information about codon usage
    print('CAI:')
    print(f'Wild type human: {calculate_cai(orf2, cu_table_human)}')
    print(f'Wild type E coli: {calculate_cai(orf2, cu_table_ecoli)}')
    print(f'Optimized E coli: {calculate_cai(best_seq, cu_table_ecoli)}')

    print('Delta MinMax fullseq:')
    print(f'Wild type E coli: {delta_min_max(minmax_native_human, minmax_native_ecoli)}')
    print(f'Optimized type E coli: {delta_min_max(minmax_native_human,min_max_full_best_seq)}')

    print('Delta MinMax harmo:')
    print(f'Wild type E coli: {delta_min_max(minmax_native_human[48:], minmax_native_ecoli[48:])}')
    print(f'Optimized type E coli: {delta_min_max(minmax_native_human[48:],min_max_full_best_seq[48:])}')

    print('Correlation MinMax fullseq:')
    print(f'Wild type E coli: {correlation_min_max(minmax_native_human, minmax_native_ecoli)}')
    print(f'Optimized type E coli: {correlation_min_max(minmax_native_human,min_max_full_best_seq)}')

    print('Correlation MinMax harmo:')
    print(f'Wild type E coli: {correlation_min_max(minmax_native_human[48:], minmax_native_ecoli[48:])}')
    print(f'Optimized type E coli: {correlation_min_max(minmax_native_human[48:],min_max_full_best_seq[48:])}')


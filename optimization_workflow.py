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
import sys
from os.path import isfile


def check_args(args):
    if len(args) != 7:
        raise TypeError("The algorithm expects 6 entries")
    if not (isfile(args[1]) and isfile(args[2]) and isfile(args[3])):
        raise TypeError("The three first parameters must be  valid filepath")
    if not all(x in ['A', 'T', 'G', 'C'] for x in args[4]):
        raise TypeError("The SD sequence must only contain A,T,G,C letters")
    if not all(x in ['A', 'T', 'G', 'C'] for x in args[5]):
        raise TypeError("The spacer sequence must only contain A,T,G,C letters")


if __name__ == '__main__':
    args = sys.argv

    check_args(args)

    CODON_WINDOW_SIZE = 10
    CODON_SHIFT = 1

    # loading codon usage tables
    cu_table_origin = CodonUsage(args[1], 'cub')       # the origin CU table
    cu_table_dest = CodonUsage(args[2], 'cub')                      # the destination CU table

    # get proteins (encoded with DNA string)
    proteins = import_proteins(args[3])
    prot = proteins[0]

    # split the seq in two parts
    # the first 48 nucleotides to lower the codon usage frequency
    prot_start_low_codons_freq = replace_codons_with_rare(prot.get_dna_seq()[:48], cu_table_dest)

    print("Sequence 16 first codons freq lowered")

    # the rest to harmonize with dest cu_table
    prot_rest = Protein(prot.get_name(), ''.join(prot.get_codon_list()[16:]), 'dna')
    prot_rest_harmonized = charming(prot_rest, cu_table_origin, cu_table_dest, 'MM', 1, CODON_WINDOW_SIZE)

    print("Sequence harmonized")

    # recombine in one single protein
    prot_full_harmonized = Protein(
        prot.get_name(),
        prot_start_low_codons_freq + prot_rest_harmonized[0].get_dna_seq(),
        'dna'
    )

    # remove alternative starts
    alt_starts = find_alternative_starts(prot_full_harmonized.get_dna_seq())
    new_seq = replace_codons(prot_full_harmonized.get_dna_seq(), alt_starts)

    print("Alternative start removed")

    prot_alt_start_removed = Protein(
        prot.get_name(),
        new_seq,
        'dna'
    )

    struct_reduction_strategy = SimulatedAnnealing(max_iteration=1000)
    struct_reduction = StructReduction(struct_reduction_strategy, cu_table_dest)

    sd = args[4]
    spacer = args[5]

    best_spacer, best_seq, best_score, iter_best = struct_reduction.reduce_structure(prot_alt_start_removed.get_dna_seq(), sd, spacer, 0.4)
    print(f'Best score {best_score}')
    print(f'found after {iter_best} iterations')
    structure, mfe = RNA.fold(sd + best_spacer + best_seq)
    print('Best sequence:')
    print(sd + best_spacer + best_seq)
    print('Best structure:')
    print(structure)

    outname = args[6]

    # Génération de la notation PostScript
    ps_file = f"{outname}.ps"
    RNA.PS_rna_plot(sd + best_spacer + best_seq, structure, ps_file)

    # Convertir la notation PostScript en image PNG
    img_file = f"rna_folding_optimized_{outname}.png"

    img = Image.open(ps_file)
    img.save(img_file, "png")

    best_seq = Protein(prot.get_name(), best_seq, 'dna')

    # calculate minmax
    minmax_native_origin = calculate_min_max(prot, CODON_WINDOW_SIZE, cu_table_origin, CODON_SHIFT)
    minmax_native_destination = calculate_min_max(prot, CODON_WINDOW_SIZE, cu_table_dest, CODON_SHIFT)
    min_max_full_best_seq = calculate_min_max(best_seq, CODON_WINDOW_SIZE, cu_table_dest, CODON_SHIFT)

    # output the result in a chart
    df = pd.DataFrame(zip(list(minmax_native_origin), list(minmax_native_destination), list(min_max_full_best_seq), ),
                      columns=['wt_origin_host', 'wt_destination_host', 'optimized_destination_host'])
    # compute the ratio
    df['reference'] = df['wt_origin_host'] - df['wt_origin_host']
    df['ratio_wild_type_destination'] = df['wt_destination_host'] - df['wt_origin_host']
    df['ratio_optimized_destination'] = df['optimized_destination_host'] - df['wt_origin_host']

    df.drop(['wt_origin_host', 'wt_destination_host', 'optimized_destination_host'], inplace=True, axis=1)

    fig = px.line(df)
    fig.update_yaxes(range=[-100, 100])
    fig.write_image(f"fig_{outname}_full_workflow.jpeg")

    # compute information about codon usage
    print('CAI:')
    print(f'Wild type origin: {calculate_cai(prot, cu_table_origin)}')
    print(f'Wild type Destination: {calculate_cai(prot, cu_table_dest)}')
    print(f'Optimized Destination: {calculate_cai(best_seq, cu_table_dest)}')

    print('Delta MinMax fullseq:')
    print(f'Wild type Destination: {delta_min_max(minmax_native_origin, minmax_native_destination)}')
    print(f'Optimized type Destination: {delta_min_max(minmax_native_origin, min_max_full_best_seq)}')

    print('Delta MinMax harmo:')
    print(f'Wild type Destination: {delta_min_max(minmax_native_origin[48:], minmax_native_destination[48:])}')
    print(f'Optimized type E Destination: {delta_min_max(minmax_native_origin[48:], min_max_full_best_seq[48:])}')

    print('Correlation MinMax fullseq:')
    print(f'Wild type Destination: {correlation_min_max(minmax_native_origin, minmax_native_destination)}')
    print(f'Optimized type Destination: {correlation_min_max(minmax_native_origin, min_max_full_best_seq)}')

    print('Correlation MinMax harmo:')
    print(f'Wild type Destination: {correlation_min_max(minmax_native_origin[48:], minmax_native_destination[48:])}')
    print(f'Optimized type Destination: {correlation_min_max(minmax_native_origin[48:], min_max_full_best_seq[48:])}')


import time
import RNA
from utils.import_proteins_file import import_proteins
from codon_usage.codon_usage import CodonUsage
from utils.protein import Protein
from structure_reduction.simulated_annealing import SimulatedAnnealing
from structure_reduction.replace_matching import ReplaceMatchingCodonsStrategy
from structure_reduction.all_possible_rna import TestAllRnasStrategy
from structure_reduction.struct_reduction import StructReduction


if __name__ == '__main__':
    # get hepatitis_e protein dna sequence
    proteins = import_proteins(file_path='test_files/hepatitis_e_virus.fasta')
    hepatitis_e = proteins[0]

    # get ecoli cu_table
    cu_table_ecoli = CodonUsage('test_files/CUB_Kazusa_Ecoli_K12.txt', 'cub')  # the destination CU table

    # create the different optimization strategies
    strategies = []
    strategies.append(TestAllRnasStrategy())
    strategies.append(SimulatedAnnealing(max_iteration=500))
    strategies.append(ReplaceMatchingCodonsStrategy(False))
    strategies.append(ReplaceMatchingCodonsStrategy(True))

    # Use of a classical SD sequence and a short spacer for the example
    sd = 'AGGAGG'
    spacer = 'A'

    # perform optimization for each strategies and generate small report in the terminal
    for strategy in strategies:
        struct_reduction = StructReduction(strategy, cu_table_ecoli)
        best_spacer, best_seq, best_score, number_iteration_best = struct_reduction.reduce_structure(hepatitis_e.get_dna_seq(), sd, spacer)
        print(f"Best spacer is {best_spacer}")
        print(f"Best score is {best_score}")
        print(f"Best score found after {number_iteration_best}")
        print("Best sequence:")
        print(sd + spacer + best_seq)
        structure, mfe = RNA.fold(sd + best_spacer + best_seq)
        print(structure)



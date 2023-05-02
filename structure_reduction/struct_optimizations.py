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
    proteins = import_proteins(file_path='test/test_files/hepatitis_e_virus.fasta')
    hepatitis_e = proteins[0]
    cu_table_scer = CodonUsage('test/test_files/CUB_Kazusa_Ecoli_K12.txt', 'cub')  # the destination CU table

    struct_reduction_strategy = SimulatedAnnealing(max_iteration=50)
    #struct_reduction_strategy = TestAllRnasStrategy()
    #struct_reduction_strategy = ReplaceMatchingCodonsStrategy(True)
    struct_reduction = StructReduction(struct_reduction_strategy, cu_table_scer)

    sd = 'AGGAGG'
    spacer = 'A'

    start = time.time()
    best_spacer, best_seq, best_score = struct_reduction.reduce_structure(hepatitis_e.get_dna_seq(), sd, spacer)
    print(best_score)
    structure, mfe = RNA.fold(sd + best_spacer + best_seq)
    print(sd + best_spacer + best_seq)
    print(structure)
    print(time.time() - start)
    print("AA sequence maintained:")
    print(hepatitis_e.get_aa_seq() == Protein("hepatitis_e", best_seq, "dna").get_aa_seq())
    print(hepatitis_e.get_aa_seq())
    print(Protein("hepatitis_e", best_seq, "dna").get_aa_seq())
    print("Nucleotides sequences conserved after 10 first codons")
    print(hepatitis_e.get_dna_seq()[len(sd+spacer)+30:297] == Protein("hepatitis_e",best_seq,"dna").get_dna_seq()[len(sd+spacer)+30:])


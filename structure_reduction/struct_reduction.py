from abc import ABC, abstractmethod
import RNA
from codon_usage.codon_usage import CodonUsage


def evaluate_structure(sequence, n_nucleotides):
    # Use RNA fold to predict secondary structure
    structure, mfe = RNA.fold(sequence)

    # Evaluate the number of nucleotides involved in a secondary structure for N_nucleotides length
    num_involved = 0
    for i in range(n_nucleotides):
        if structure[i] == '(':
            num_involved += 1
    return num_involved, structure


class StructReductionStrategy(ABC):
    @abstractmethod
    def run(self, sequence, sd, spacer):
        pass


class StructReduction():
    def __init__(self, struct_reduction_strategy: StructReductionStrategy, cu_table: CodonUsage):
        self._struct_reduction_strategy = struct_reduction_strategy
        self._cu_table = cu_table

    def reduce_structure(self, sequence, sd, spacer, threshold):
        return self._struct_reduction_strategy.run(sequence, sd, spacer, self._cu_table, threshold)
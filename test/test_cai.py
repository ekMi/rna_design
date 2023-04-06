from utils.protein import Protein
from codon_usage.codon_usage import CodonUsage
from codon_usage.cai import calculate_cai

if __name__ == '__main__':
    cu_table = CodonUsage('test_files/VeryHighLevelExpressionEColi.txt')
    protein_file = open('test_files/test_proteins.txt', 'r')

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

    print(len(proteins))

    for protein in proteins:
        print(protein.get_name())
        print(round(calculate_cai(protein, cu_table),3))


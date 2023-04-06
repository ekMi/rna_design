from utils.protein import Protein


def import_proteins(file_path='test_files/test_proteins.txt'):
    protein_file = open(file_path, 'r')

    lines = protein_file.readlines()

    proteins = []
    for count, line in enumerate(lines):
        if ">" in line:
            if count != 0:
                proteins.append(Protein(protein_name, sequence, 'dna'))
            sequence = ""
            protein_name = line.strip().removeprefix(">")
        else:
            sequence += line.strip()
    proteins.append(Protein(protein_name, sequence, 'dna'))

    return proteins

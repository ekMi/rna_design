import sys
from protein import Protein

mapDict = {'TCA': 'S', 'AAT': 'N', 'TGG': 'W', 'GAT': 'D', 'GAA': 'E', 'TTC': 'F', 'CCG': 'P',
           'ACT': 'T', 'GGG': 'G', 'ACG': 'T', 'AGA': 'R', 'TTG': 'L', 'GTC': 'V', 'GCA': 'A',
           'TGA': '*', 'CGT': 'R', 'CAC': 'H', 'CTC': 'L', 'CGA': 'R', 'GCT': 'A', 'ATC': 'I',
           'ATA': 'I', 'TTT': 'F', 'TAA': '*', 'GTG': 'V', 'GCC': 'A', 'GAG': 'E', 'CAT': 'H',
           'AAG': 'K', 'AAA': 'K', 'GCG': 'A', 'TCC': 'S', 'GGC': 'G', 'TCT': 'S', 'CCT': 'P',
           'GTA': 'V', 'AGG': 'R', 'CCA': 'P', 'TAT': 'Y', 'ACC': 'T', 'TCG': 'S', 'ATG': 'M',
           'TTA': 'L', 'TGC': 'C', 'GTT': 'V', 'CTT': 'L', 'CAG': 'Q', 'CCC': 'P', 'ATT': 'I',
           'ACA': 'T', 'AAC': 'N', 'GGT': 'G', 'AGC': 'S', 'CGG': 'R', 'TAG': '*', 'CGC': 'R',
           'AGT': 'S', 'CTA': 'L', 'CAA': 'Q', 'CTG': 'L', 'GGA': 'G', 'TGT': 'C', 'TAC': 'Y',
           'GAC': 'D'}

aaDict = {'S': ['TCA', 'TCC', 'TCT', 'TCG', 'AGC', 'AGT'], 'N': ['AAT', 'AAC'], 'W': ['TGG'],
          'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'F': ['TTC', 'TTT'], 'P': ['CCG', 'CCT', 'CCA', 'CCC'],
          'T': ['ACT', 'ACG', 'ACC', 'ACA'], 'G': ['GGG', 'GGC', 'GGT', 'GGA'],
          'R': ['AGA', 'CGT', 'CGA', 'AGG', 'CGG', 'CGC'], 'L': ['TTG', 'CTC', 'TTA', 'CTT', 'CTA', 'CTG'],
          'V': ['GTC', 'GTG', 'GTA', 'GTT'], 'A': ['GCA', 'GCT', 'GCC', 'GCG'], '*': ['TGA', 'TAA', 'TAG'],
          'H': ['CAC', 'CAT'], 'I': ['ATC', 'ATA', 'ATT'], 'K': ['AAG', 'AAA'], 'Y': ['TAT', 'TAC'],
          'M': ['ATG'], 'C': ['TGC', 'TGT'], 'Q': ['CAG', 'CAA']}

def generateAAFreqDict(freqDict, aaDict):
    """
    Takes in a dictionary mapping codons to frequencies and amino acids to codons
    and returns a dictionary mapping amino acids to a list of theirrespective
    codon frequencies
    """
    aaFreqDict = dict()

    for aa in aaDict:
        aaFreqDict[aa] = []
        for codon in aaDict[aa]:
            aaFreqDict[aa].append(freqDict[codon])

    return aaFreqDict

def calculateMinMax(sequence, aaFreqDict, freqDict, mapDict, windowSize):
    """This function calculates %MinMax values for a given sequence, returned as a list of floats"""

    minMaxValues = []

    for i in range(int(windowSize / 2)):
        minMaxValues.append(0)

    # Using the specified sliding window size (windowSize/2 - 1 on either side of the central codon), min/max is calculated
    for i in range(len(sequence) - windowSize + 1):
        window = sequence[i:i + windowSize]  # list of the codons in the current window

        Actual = 0.0  # average of the actual codon frequencies
        Max = 0.0  # average of the min codon frequencies
        Min = 0.0  # average of the max codon frequencies
        Avg = 0.0  # average of the averages of all the frequencies associated with each amino acid

        # Sum the frequencies
        for codon in window:
            frequencies = aaFreqDict[
                mapDict[codon]]  # list of all frequencies associated with the amino acid this codon encodes
            Actual += freqDict[codon]
            Max += max(frequencies)
            Min += min(frequencies)
            Avg += sum(frequencies) / len(frequencies)

        # Divide by the window size to get the averages
        Actual = Actual / windowSize
        Max = Max / windowSize
        Min = Min / windowSize
        Avg = Avg / windowSize

        percentMax = ((Actual - Avg) / (Max - Avg)) * 100
        percentMin = ((Avg - Actual) / (Avg - Min)) * 100

        if (percentMax >= 0):
            minMaxValues.append(round(percentMax, 2))
        else:
            minMaxValues.append(round(-percentMin, 2))

    # fills in values for codons where window size makes min/max unable to be calculated
    if windowSize % 2 == 1:
        for i in range(int(windowSize / 2)):
            minMaxValues.append(0)
    else:
        for i in range(int(windowSize / 2) - 1):
            minMaxValues.append(0)
    return minMaxValues

def readCUBTable(path):
    """
    This function reads in a CUB table where the format for a line is:
    <Codon> <Frequency>
    repeated for each codon.
    """
    cubDict = dict()
    cubTable = open(path, "r")
    for line in cubTable:
        line = line.replace("U", "T")
        line = line.split()
        if len(line) == 2:
            if line[0] in mapDict:
                cubDict[line[0]] = float(line[1])

    if len(cubDict) != len(mapDict):
        print("The following codons are missing/in the incorrect format in your input CUB table:")
        for codon in mapDict:
            if codon not in cubDict:
                print(codon)
        print("Please fix these and rerun the script.")
        sys.exit()
    return cubDict

if __name__ == '__main__':
    freqDict = readCUBTable('./test_files/CUB_HiveCut_Ecoli_K12.txt')
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
    proteins.append(Protein(protein_name, sequence, 'dna'))

    # test on the third protein (rpoD) which is relatively long
    rpoD = proteins[2]
    aaFreqDict = generateAAFreqDict(freqDict, aaDict)
    print(calculateMinMax(rpoD.get_codon_list(), aaFreqDict, freqDict, mapDict, 18))
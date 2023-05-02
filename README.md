# Optimization Workflow

## Description

This program performs an optimization workflow to
- harmonize the codon usage frequency from a wild type host to a heterlogous host
- reduce the codon usage frequency at the beginning of the sequence
- remove alternative starts
- reduce the secondary structure at the beginning of the sequence using a simulated annealing strategy

## Prerequisites

- Python 3.9
- The python dependies of the requirements.txt file
- Tutorial files can be used to try the algorithm and are located in the `tutorial_files` directory:
    - codon usage table for the origin host (`cu_table_origin.txt`)
    - codon usage table for the destination host (`cu_table_dest.txt`)
    - protein sequence(s) in FASTA format (`proteins.fasta`)

## Usage

To run the program, execute the following command in the terminal:

```
python optimization_workflow.py <path_to_cu_table_origin> <path_to_cu_table_dest> <path_to_protein_fasta_file> <SD_sequence> <spacer_sequence> <output_name>
```

Where:
- `<path_to_cu_table_origin>` is the path to the codon usage table for the origin host
- `<path_to_cu_table_dest>` is the path to the codon usage table for the destination host
- `<path_to_protein_fasta_file>` is the path to the FASTA file containing the protein sequence(s)
- `<SD_sequence>` is the Shine-Dalgarno sequence to use for translation initiation
- `<spacer_sequence>` is the spacer sequence to use between the Shine-Dalgarno sequence and the start codon
- `<output_name>` is the name to use for the output files

## Output

The program generates the following output:

- `output_name.png`: a PNG image of the RNA secondary structure
- `output_name.ps`: a PostScript file of the RNA secondary structure
- a chart showing the min-max scores for the wild-type sequence (in the origin host), the wild-type sequence (in the destination host), and the optimized sequence (in the destination host)
- several outputs in the terminal like CAI of the different sequences, the correlation minmax,...
from PIL import Image
import RNA
from utils.import_proteins_file import import_proteins

# Séquence d'ARN

proteins = import_proteins("test_files/hepatitis_e_virus.fasta")
orf2 = proteins[0]

seq = orf2.get_dna_seq()

sd = 'AGGAGG'
spacer = 'AAAAAAAA'

# Calcul de la structure en notation point-parenthèse
structure = RNA.fold(sd + spacer + seq)[0]
print(structure)

# Génération de la notation PostScript
ps_file = "orf2_native.ps"
RNA.PS_rna_plot(sd + spacer + seq, structure, ps_file)

# Convertir la notation PostScript en image PNG
img_file = "rna_folding_native_orf2.png"
img = Image.open(ps_file)
img.save(img_file, "png")

# Afficher l'image
img.show()

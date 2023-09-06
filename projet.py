from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices

def blosum_score():
    # Accédez à la matrice BLOSUM62
    blosum62 = substitution_matrices.load("BLOSUM62")

    # Définissez les pénalités pour les gaps
    gap_open = -1
    gap_extend = -1

    # Utilisez pairwise2 pour effectuer l'alignement avec BLOSUM62 et les pénalités de gap
    alignments = pairwise2.align.globalds("TRAVAIL", "CARNAVAL", blosum62, gap_open, gap_extend, one_alignment_only=True)

    # Obtenez le meilleur alignement
    best_alignment = alignments[0]

    # Affichez l'alignement formaté
    print("Alignment:")
    print(format_alignment(*best_alignment))

if __name__ =="__main__":
    blosum_score()

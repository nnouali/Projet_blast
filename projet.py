import numpy as np
from Bio.Align import substitution_matrices


def blosum_score():
    # Accédez à la matrice BLOSUM62 à l'aide de Bio.Align.substitution_matrices
    blosum62 = substitution_matrices.load("BLOSUM62")

    # Exemple d'accès aux scores de substitution
    score_LD = blosum62["L"]["D"]  # Score de substitution entre A et R
    print(f"Score de substitution entre L et D: {score_LD}")


if __name__ =="__main__":
    blosum_score()

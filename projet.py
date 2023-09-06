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
'''
def calculate_score(query, database_sequence):
    # Calculez le score de correspondance entre la séquence de requête et la séquence de la base de données.
    score = 0
    for i in range(len(query)):
        if query[i] == database_sequence[i]:
            score += 1
    return score

def blast(query, database, threshold):
    results = []
    for db_sequence in database:
        for i in range(len(db_sequence) - len(query) + 1):
            subsequence = db_sequence[i:i+len(query)]
            score = calculate_score(query, subsequence)
            if score >= threshold:
                results.append((db_sequence, i, score))
    return results

# Exemple d'utilisation
query_sequence = "AGTCA"
database_sequences = ["AGTCCATGC", "AGTCACTGA", "TGGGAGTCA"]
threshold_score = 17

results = blast(query_sequence, database_sequences, threshold_score)

for result in results:
    print(f"Alignement trouvé dans la séquence : {result[0]} à la position {result[1]} avec un score de {result[2]}")
'''
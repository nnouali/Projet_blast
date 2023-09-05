import numpy as np

# matrice de score avec les bons para
def score_matrix (rows, colums) :
    matrix = [] 
    for i in range (rows):
        row = []
        for j in range (colums):
            score = 0
            row.apped(score)
        matrix.append()
    return matrix
#ou 
matrix_similarite = np.array([
    [0, 0, 30],
    [50, 5, 25],
    [8, 12, 22]
])

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

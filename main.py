from Bio.Align import substitution_matrices
from Bio import SeqIO
import sqlite3

from database import split_sequence_into_words, read_fasta_file

from projet import FindWords, CalculateScore, FindNeighbors, CompareNeighborsWithDatabase, blosum62Map, ExtendWordInDatabase

blosum62 = substitution_matrices.load("BLOSUM62")

if __name__ == "__main__":
    # Lire le fichier FASTA
    fasta_file = "./test.fasta"  # Remplacez par le chemin de votre fichier FASTA
    sequence = str(SeqIO.read(fasta_file, "fasta").seq)

    word_length = 3
    threshold = 11

    # Trouver les mots de la séquence
    words = FindWords(sequence, word_length)

    # Trouver les voisins avec leurs scores
    neighbors = FindNeighbors(words, word_length, threshold)

    # Comparer les voisins avec la base de données
    database_path = 'sequence_database.db'  # Remplacez par le chemin de votre base de données
    matching_neighbors = CompareNeighborsWithDatabase(neighbors, database_path)

    # Ouvrir le fichier de sortie
    output_file = "mots_voisins_scores.txt"
    with open(output_file, 'w') as file:
        for neighbor in matching_neighbors:
            sequence_id, neighbor_word, original_word, score, position = neighbor
            file.write(f"ID de séquence : {sequence_id}, Mot Voisin : {neighbor_word}, Mot d'Origine : {original_word}, Score : {score}, Position : {position}\n")
    
    print("Comparaison avec la base de données terminée.")

    # Parcourir les mots correspondants et afficher l'alignement avec les scores cumulés
    for neighbor in matching_neighbors:
        sequence_id, neighbor_word, original_word, score, position = neighbor
        print(f"ID de séquence : {sequence_id}, Mot Voisin : {neighbor_word}, Mot d'Origine : {original_word}, Score : {score}, Position : {position}")

        # Utiliser la fonction pour aligner le mot dans la séquence de la base de données
        alignment, max_score = ExtendWordInDatabase(sequence, neighbor_word, position, blosum62Map)
        print(f"Alignement : {alignment}, Score maximal : {max_score}")



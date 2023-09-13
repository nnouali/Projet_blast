from Bio.Align import substitution_matrices
from Bio import SeqIO
import sqlite3
import difflib

from database import (
    split_sequence_into_words,
    read_fasta_file,
    conn,
    cur,
    insert_sequence_if_not_exists,
)
from projet import (
    FindWords,
    FindNeighbors,
    CompareNeighborsWithDatabase,
    ExtendSequenceWithNeighbor,
    blosum62Map,
)

blosum62 = substitution_matrices.load("BLOSUM62")

if __name__ == "__main__":
    # Initialisation des variables pour stocker le meilleur score et l'alignement correspondant
    best_score = None
    best_alignment_info = None
    # Lire le fichier FASTA
    read_fasta_file("Liste_Fasta2.txt")

    fasta_file = "./test.fasta"  # Remplacez par le chemin de votre fichier FASTA
    sequence = str(SeqIO.read(fasta_file, "fasta").seq)
    query_sequence = str(SeqIO.read(fasta_file, "fasta").seq)
    # Création d'une connexion à la base de données
    conn = sqlite3.connect('sequence_database.db')
    cur = conn.cursor()

    # Ajoutez la séquence à la base de données si elle n'existe pas déjà
    insert_sequence_if_not_exists("sequence_id", sequence)  # Remplacez "sequence_id_unique" par un identifiant unique pour la séquence

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

    # Étendre la séquence avec chaque voisin
    extended_sequences = []
    for neighbor in matching_neighbors:
        extended_sequence = ExtendSequenceWithNeighbor(query_sequence, neighbor, word_length)
        extended_sequences.append((neighbor, extended_sequence))

    # Ouvrir le fichier de sortie
    output_file = "seq_voisins_scores.txt"
    with open(output_file, 'w') as file:
        for neighbor, extended_sequence in extended_sequences:
            sequence_id, neighbor_word, original_word, score, position = neighbor
            file.write(f"ID de séquence : {sequence_id}, Mot Voisin : {neighbor_word}, Mot d'Origine : {original_word}, Score : {score}, Position : {position}, Mot Étendu : {extended_sequence}\n")

    print("Comparaison avec la base de données terminée.")

    # Créez un dictionnaire pour stocker les mots étendus avec leurs alignements et scores
    extended_word_data = {}

    with open("seq_voisins_scores.txt", 'r') as file:
        for line in file:
            parts = line.strip().split(',')
            sequence_id = parts[0].split(':')[1].strip()
            db_sequence_id = parts[0].split(':')[1].strip()
            neighbor_word = parts[1].split(':')[1].strip()
            extended_word = parts[-1].split(':')[1].strip()

            # Récupérez la séquence de la base de données correspondante
            cur.execute('SELECT sequence FROM sequences WHERE sequence_id = ?', (sequence_id,))
            db_sequence = cur.fetchone()[0]

            # Calculez le score par position et créez une ligne d'alignement avec des scores
            alignment = []
            scores = []
            for i in range(min(len(extended_word), len(db_sequence))):
                if extended_word[i] == db_sequence[i]:
                    alignment.append('|')  # Correspondance
                else:
                    alignment.append(' ')
                # Calculez le score et stockez-le
                score = blosum62Map[(extended_word[i], db_sequence[i])]
                scores.append(score)

            # Créez une ligne d'alignement avec des traits
            alignment = ''.join('|' if extended_word[i] == db_sequence[i] else ' ' for i in range(min(len(extended_word), len(db_sequence))))

            # Stockez les données
            if (sequence_id, neighbor_word) not in extended_word_data:
                extended_word_data[(sequence_id, neighbor_word)] = {
                    'extended_word': extended_word,
                    'alignment': alignment,
                    'db_sequence': db_sequence,
                    'scores': scores
                }

    # Affichez l'alignement des mots étendus avec les scores à chaque position
    for (sequence_id, neighbor_word), data in extended_word_data.items():
        extended_word = data['extended_word']
        alignment = data['alignment']
        db_sequence = data['db_sequence']
        scores = data['scores']
        # Récupérez le mot d'origine à partir des données de chaque voisin
        matching_neighbor = [neighbor for neighbor in matching_neighbors if neighbor[1] == neighbor_word][0]
        original_word = matching_neighbor[2]

        print()
        print(f"Request sequence ID: {sequence_id}, Mot Voisin : {neighbor_word}")
        print("Request Word      : ", original_word)
        print(f"Extended Word : {extended_word}")
        print(f" Database sequence : {db_sequence}")
        print("Alignement :")
        print("     Extended Word     : ", extended_word)
        print("                       : ", alignment)
        print("     Database sequence : ", db_sequence)
        print("     Scores            : ", ' '.join([f"{score:.1f}" for score in scores]))

        # Calculer les scores cumulés
        cumulative_scores = []
        cumulative_score = 0.0
        for score in scores:
            cumulative_score += score
            cumulative_scores.append(cumulative_score)

        print("     Scores Cumulés    : ", ' '.join([f"{score:.1f}" for score in cumulative_scores]))

        # Trouver l'indice k où le mot voisin correspondant à db_sequence commence
        k = db_sequence.find(neighbor_word)

        if k != -1:
            # Calculer les scores cumulés à partir de la position k vers la droite
            cumulative_scores_d = []
            cumulative_score_d = 0.0
            for l in range(k, len(scores)):
                cumulative_score_d += scores[l]
                cumulative_scores_d.append(cumulative_score_d)
            
            print("     Cumulative score RIGHT : ", ' '.join([f"{score:.1f}" for score in cumulative_scores_d]))
        else:
            print("     Cumulative score RIGHT DOES NOT EXIST")

        # Trouver l'indice i où le mot voisin correspondant à db_sequence commence
        i = db_sequence.find(neighbor_word)

        if i != -1:
            # Commencer à cumuler à partir de la fin du mot voisin
            start_position = i + len(neighbor_word) - 1  # Position du dernier acide aminé du mot voisin dans db_sequence
            
            # Calculer les scores cumulés à partir de la position start_position vers la gauche
            cumulative_scores_g = []
            cumulative_score_g = 0.0
            for j in range(start_position, -1, -1):  # Commencez à partir de start_position jusqu'au début (de droite à gauche)
                cumulative_score_g += scores[j]
                cumulative_scores_g.insert(0, cumulative_score_g)  # Insérer au début de la liste
                
            print("     Cumulative score LEFT : ", ' '.join([f"{score:.1f}" for score in cumulative_scores_g]))
        else:
            print("     Cumulative score LEFT DOES NOT EXIST")

        # Calculer le maximum du Scores Cumulés à partir du Mot Voisin
        max_cumulative_score_g = max(cumulative_scores_g)
        max_cumulative_score_d = max(cumulative_scores_d)

        # Trouver les indices où les maximums des Scores Cumulés à partir du Mot Voisin sont atteints
        max_cumulative_score_index_g = cumulative_scores_g.index(max_cumulative_score_g)
        max_cumulative_score_index_d = cumulative_scores_d.index(max_cumulative_score_d)

        # Trouver l'indice k où le mot voisin correspondant à db_sequence commence
        k = db_sequence.find(neighbor_word)

        # Extraire la séquence à gauche et à droite du maximum
        sequence_left = extended_word[max_cumulative_score_index_g:k + len(neighbor_word)-3]
      # Trouver l'indice k où le mot voisin correspondant à db_sequence commence
        n = db_sequence.find(neighbor_word)

        # Extraire la séquence de droite à partir du neighbor_word jusqu'au max_cumulative_score
        sequence_right = extended_word[n:n + max_cumulative_score_index_d + 1]


        # Combinez les deux parties pour obtenir l'alignement complet
        alignment_complete = sequence_left + sequence_right
        
        print("\nFinal Alignement :")
        # Afficher l'alignement complet
        print("     Request sequence  :", alignment_complete)

        # Extraire la séquence à gauche et à droite de la db_sequence correspondante
        db_sequence_left = db_sequence[max_cumulative_score_index_g:k + len(neighbor_word)-3]

        # Extraire la séquence à droite de la db_sequence correspondante
        db_sequence_right = db_sequence[n:n + max_cumulative_score_index_d + 1]

        alignment_complete_db = db_sequence_left + db_sequence_right

        # Afficher l'alignement complet
        print("                      :", alignment)
        print("     Sequence database :", alignment_complete_db)

        # Calculer le score entre alignment_complete et alignment_complete_db
        new_score = sum(blosum62Map[pair] for pair in zip(alignment_complete, alignment_complete_db))

        # Afficher le nouveau score
        print("\nFinal score:", new_score)

        # Comparer le nouveau score avec le meilleur score actuel
        if best_score is None or new_score > best_score:
            best_score = new_score

            # Construire l'information d'alignement complète
            best_alignment_info = [
                "*************************************************************************",
                "\n* Best Alignement : ",
                "*",
                "*   Request sequence     : " + alignment_complete,
                "*                        :" + alignment,
                "*   Sequence database    : " + alignment_complete_db,
                "*   Request sequence ID  :" + sequence_id,
                "*   Database sequence ID :" + db_sequence_id,               
                "**************************************************************************",
            ]
            
    # Afficher le meilleur score et l'alignement correspondant
    if best_score is not None:
        print("\nBEST SCORE :", best_score)
        print("\n".join(best_alignment_info))
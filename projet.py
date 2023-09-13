from Bio.Align import substitution_matrices
import sqlite3

# Charger la matrice de substitution BLOSUM62
blosum62 = substitution_matrices.load("BLOSUM62")



# Charger la matrice de substitution BLOSUM62
blosum62 = substitution_matrices.load("BLOSUM62")

aminoacids = ['A', 'C', 'D', 'E',
              'F', 'G', 'H', 'I',
              'K', 'L', 'M', 'N',
              'P', 'Q', 'R', 'S',
              'T', 'V', 'W', 'Y']

blosum62Map = {}
for i in aminoacids:
    for j in aminoacids:
        blosum62Map[(i, j)] = blosum62[i][j]


def FindWords(inputSequence, wordLength):
    wordsList = []
    for i in range(0, len(inputSequence) - wordLength + 1):
        outputWord = inputSequence[i: i + wordLength]
        wordsList.append((outputWord, i))
    return wordsList

def CalculateScore(word1, word2, blosum62Map):
    score = sum(blosum62Map[(word1[i], word2[i])] for i in range(len(word1)))
    return score

def FindNeighbors(words, wordLength, threshold):
    def generate_neighbors(word, depth, wordIndex):
        if depth == wordLength:
            finalScore = CalculateScore(word, words[wordIndex][0], blosum62Map)
            if finalScore >= threshold:
                neighborsSet.add((word, words[wordIndex][0], finalScore, words[wordIndex][1]))
            return
        for i in range(wordLength):
            for aminoAcid in aminoacids:
                if aminoAcid != word[i]:
                    neighborWord = word[:i] + aminoAcid + word[i+1:]
                    generate_neighbors(neighborWord, depth + 1, wordIndex)

    neighborsSet = set()
    for (word, wordIndex) in words:
        generate_neighbors(word, 0, wordIndex)

    neighborsList = sorted(neighborsSet, key=lambda x: x[2], reverse=True)
    return neighborsList

def SearchWordsInDatabase(sequence, word_length, threshold, conn, cur):
    words = set()
    for i in range(0, len(sequence) - word_length + 1):
        word = sequence[i:i+word_length]
        words.add(word)
    
    results = {}
    
    for word in words:
        cur.execute('''
            SELECT sequence_id, position
            FROM sequence_words
            WHERE word = ?
        ''', (word,))
        
        matching_words = cur.fetchall()
        
        for sequence_id, position in matching_words:
            results[(sequence_id, word)] = position
    
    return results

def CompareNeighborsWithDatabase(neighbors, database_path):
    # Charger la base de données SQLite
    conn = sqlite3.connect(database_path)
    cur = conn.cursor()

    # Lire les mots de la base de données avec les identifiants de séquence correspondants
    cur.execute('SELECT sequence_id, word, position FROM sequence_words')
    db_words = cur.fetchall()

    # Créer un ensemble pour les mots de la base de données
    db_word_set = {(sequence_id, word, position) for sequence_id, word, position in db_words}

    # Liste pour stocker les correspondances
    matching_neighbors = []

    for neighbor in neighbors:
        neighbor_word, original_word, score, position = neighbor
        # Vérifier si le voisin est dans la base de données
        for sequence_id, word, db_position in db_word_set:
            if word == neighbor_word and position == db_position:
                matching_neighbors.append((sequence_id, neighbor_word, original_word, score, position))

    # Fermer la connexion à la base de données
    conn.close()

    return matching_neighbors

def ExtendSequenceWithNeighbor(sequence, neighbor, word_length):
    sequence_id, neighbor_word, original_word, score, position = neighbor
    extended_sequence = sequence[:position] + neighbor_word + sequence[position + word_length:]
    return extended_sequence
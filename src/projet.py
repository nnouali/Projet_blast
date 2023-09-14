from Bio.Align import substitution_matrices
import sqlite3

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
    """
    Trouve tous les mots d'une certaine longueur dans une séquence d'entrée.

    Args:
        inputSequence (str): La séquence d'entrée.
        wordLength (int): La longueur des mots à rechercher.

    Returns:
        list: Une liste de tuples contenant les mots trouvés et leur position dans la séquence.
    """
    wordsList = []
    for i in range(0, len(inputSequence) - wordLength + 1):
        outputWord = inputSequence[i: i + wordLength]
        wordsList.append((outputWord, i))
    return wordsList


def CalculateScore(word1, word2, blosum62Map):
    """
    Calcule le score entre deux mots en utilisant la matrice de substitution BLOSUM62.

    Args:
        word1 (str): Le premier mot.
        word2 (str): Le deuxième mot.
        blosum62Map (dict): La carte des scores de la matrice de substitution BLOSUM62.

    Returns:
        int: Le score calculé.
    """
    score = sum(blosum62Map[(word1[i], word2[i])] for i in range(len(word1)))
    return score


def FindNeighbors(words, wordLength, threshold):
    """
    Trouve les voisins de mots similaires dans une liste de mots.

    Args:
        words (list): Une liste de tuples contenant des mots et leur position dans la séquence.
        wordLength (int): La longueur des mots.
        threshold (int): Le seuil de score pour considérer un voisin comme similaire.

    Returns:
        list: Une liste de voisins similaires avec leurs informations.
    """
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
    """
    Recherche les mots dans une base de données en utilisant une séquence donnée.

    Args:
        sequence (str): La séquence à rechercher.
        word_length (int): La longueur des mots.
        threshold (int): Le seuil de score pour considérer un mot comme similaire.
        conn: Connexion à la base de données.
        cur: Curseur de la base de données.

    Returns:
        dict: Un dictionnaire des résultats avec les séquences correspondantes et leur position.
    """
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
    """
    Compare les voisins avec les mots dans une base de données.

    Args:
        neighbors (list): Une liste de voisins similaires avec leurs informations.
        database_path (str): Le chemin de la base de données.

    Returns:
        list: Une liste de voisins similaires qui existent dans la base de données.
    """
    # Charge la base de données SQLite
    conn = sqlite3.connect(database_path)
    cur = conn.cursor()

    # Lit les mots de la base de données avec les identifiants de séquence correspondants
    cur.execute('SELECT sequence_id, word, position FROM sequence_words')
    db_words = cur.fetchall()

    # Crée un ensemble pour les mots de la base de données
    db_word_set = {(sequence_id, word, position) for sequence_id, word, position in db_words}

    # Liste pour stocker les correspondances
    matching_neighbors = []

    for neighbor in neighbors:
        neighbor_word, original_word, score, position = neighbor
        # Vérifie si le voisin est dans la base de données
        for sequence_id, word, db_position in db_word_set:
            if word == neighbor_word and position == db_position:
                matching_neighbors.append((sequence_id, neighbor_word, original_word, score, position))

    # Fermer la connexion à la base de données
    conn.close()

    return matching_neighbors


def ExtendSequenceWithNeighbor(sequence, neighbor, word_length):
    """
    Étend la séquence avec un voisin similaire.

    Args:
        sequence (str): La séquence d'origine.
        neighbor (tuple): Informations sur le voisin similaire.
        word_length (int): La longueur des mots.

    Returns:
        str: La séquence étendue avec le voisin similaire.
    """
    sequence_id, neighbor_word, original_word, score, position = neighbor
    extended_sequence = sequence[:position] + neighbor_word + sequence[position + word_length:]
    return extended_sequence

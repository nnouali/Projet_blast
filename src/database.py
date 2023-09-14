import sqlite3

# Connexion à la base de données
conn = sqlite3.connect('sequence_database.db')
cur = conn.cursor()

# Créez la table "sequence_words" si elle n'existe pas déjà
cur.execute('''
    CREATE TABLE IF NOT EXISTS sequence_words (
        sequence_id TEXT,
        word TEXT,
        position INTEGER
    )
''')

# Créez la table "sequences" si elle n'existe pas déjà
cur.execute('''
    CREATE TABLE IF NOT EXISTS sequences (
        sequence_id TEXT PRIMARY KEY,
        sequence TEXT
    )
''')

conn.commit()

# Créez un dictionnaire pour stocker les mots déjà insérés
inserted_words = {}

def split_sequence_into_words(sequence, sequence_id):
    """
    Split a sequence into words of a given length and insert them into the database.

    Args:
        sequence (str): The input sequence to split.
        sequence_id (str): The identifier for the sequence.

    Returns:
        None
    """
    word_length = 3
    for i in range(0, len(sequence) - word_length + 1):
        word = sequence[i:i+word_length]
        
        # Vérifiez si le mot n'existe pas déjà dans le dictionnaire
        if word not in inserted_words:
            # Insérez le mot dans la table sequence_words
            cur.execute('''
                INSERT INTO sequence_words (sequence_id, word, position)
                VALUES (?, ?, ?)
            ''', (sequence_id, word, i))
            conn.commit()
            
            # Ajoutez le mot au dictionnaire pour éviter les doublons
            inserted_words[word] = None


def insert_sequence_if_not_exists(sequence_id, sequence):
    """
    Insert a sequence into the "sequences" table if it doesn't already exist.

    Args:
        sequence_id (str): The identifier for the sequence.
        sequence (str): The sequence to insert.

    Returns:
        None
    """
    cur.execute('SELECT sequence_id FROM sequences WHERE sequence_id = ?', (sequence_id,))
    result = cur.fetchone()
    if result is None:
        cur.execute('INSERT INTO sequences (sequence_id, sequence) VALUES (?, ?)', (sequence_id, sequence))
        conn.commit()

def read_fasta_file(file_path):
    """
    Read sequences from a FASTA file and insert them into the database.

    Args:
        file_path (str): The path to the FASTA file.

    Returns:
        None
    """
    with open(file_path, 'r') as fasta_file:
        fasta_id = ""
        sequence = ""
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    insert_sequence_if_not_exists(fasta_id, sequence)
                fasta_id = line[1:]
                sequence = ""
            else:
                sequence += line
        if sequence:
            insert_sequence_if_not_exists(fasta_id, sequence)
            
def read_fasta_file_request(file_path):
    """
    Read sequences from a FASTA file, split them into words, and insert them into the database.

    Args:
        file_path (str): The path to the FASTA file.

    Returns:
        None
    """
    with open(file_path, 'r') as fasta_file:
        fasta_id = ""
        sequence = ""
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    split_sequence_into_words(sequence, fasta_id)
                fasta_id = line[1:]
                sequence = ""
            else:
                sequence += line
        if sequence:
            split_sequence_into_words(sequence, fasta_id)

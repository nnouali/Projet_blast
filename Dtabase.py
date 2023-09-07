import sqlite3

# Création d'une connexion à la base de données
conn = sqlite3.connect('sequence_database.db')

# Curseur pour interagir avec la base de données
cur = conn.cursor()
""" TODO REPRENDRE COMMENT CHARGER LA DB
# Table qui stockera les mots de trois lettres et leurs positions
cur.execute('''
    CREATE TABLE sequence_words (
        sequence_id TEXT,
        word TEXT,
        position INTEGER
    )
''')
"""
# Valider la transaction
conn.commit()

def split_sequence_into_words(sequence, sequence_id):
    """Découpe une séquence en mots de trois lettres et les enregistre avec leurs positions.

    Args:
        sequence : La séquence.
        sequence_id : L'identifiant de la séquence.
    """
    words = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    for i, word in enumerate(words):
        # Insérer le mot et sa position dans la table de la base de données
        cur.execute('''
            INSERT INTO sequence_words (sequence_id, word, position)
            VALUES (?, ?, ?)
        ''', (sequence_id, word, i))
    # Valider la transaction
    conn.commit()

def read_fasta_file(file_path):
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

if __name__ == "__main__":
    # Lisez le fichier FASTA et effectuez les opérations de découpage en mots de trois lettres
    read_fasta_file("Liste_Fasta.txt")

    conn.close()

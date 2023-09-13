import sqlite3

conn = sqlite3.connect('sequence_database.db')
cur = conn.cursor()

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

def split_sequence_into_words(sequence, sequence_id):
    word_length = 3
    for i in range(0, len(sequence) - word_length + 1):
        word = sequence[i:i+word_length]
        cur.execute('''
            INSERT INTO sequence_words (sequence_id, word, position)
            VALUES (?, ?, ?)
        ''', (sequence_id, word, i))
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

from http.client import LineTooLong
import sqlite3

"""
# Création d'une connexion à la base de données 
conn = sqlite3.connect('ma_base_de_donnees.db')

# Curseur pour interagir avec la base de données
cursor = conn.cursor()

# Table qui stocke les séquences
cursor.execute('''
    CREATE TABLE sequences (
        id INTEGER PRIMARY KEY,
        nom TEXT,
        description TEXT,
        sequence TEXT
    )
''')"""

def read_fasta(fasta_text):
    sequences = []
    fasta_id = ""
    description = ""
    sequence = ""

    for line in fasta_text.split('\n'):
        line = line.strip()
        if line.startswith(">"):
            if sequence:
                sequences.append((fasta_id, description, sequence))
            fasta_id, description = line[1:].split(" ", 1)
            sequence = ""
        else:
            sequence += line

    if sequence:
        sequences.append((fasta_id, description, sequence))
        
    return sequences

if __name__ == "__main__":
    fasta_data = ">sp|A0A0C5B5G6|MOTSC_HUMAN Mitochondrial-derived peptide MOTS-c OS=Homo sapiens OX=9606 GN=MT-RNR1 PE=1 SV=1\nMRWQEMGYIFYPRKLR"
    sequences = read_fasta(fasta_data)

    for fasta_id, description, sequence in sequences:
        print("Nom:", fasta_id)
        print("Description:", description)
        print("Séquence:", sequence)

#sp|A0A0C5B5G6|MOTSC_HUMAN Mitochondrial-derived peptide MOTS-c OS=Homo sapiens OX=9606 GN=MT-RNR1 PE=1 SV=1
#MRWQEMGYIFYPRKLR






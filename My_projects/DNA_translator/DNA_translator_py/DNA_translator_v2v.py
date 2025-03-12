from Bio import SeqIO
from Bio.Seq import Seq
import time

# Start time
start_time = time.time()

def DNA_translator(file_FASTA):
    # Load DNA sequence from FASTA file
    with open(fasta_file, "r") as file:
        record = next(SeqIO.parse(file, "fasta"))

    # Get the sequence
    dna_seq = record.seq
    print("Loaded DNA Sequence:", dna_seq)

    # Get Reverse Complement
    reverse_complement = dna_seq.reverse_complement()
    print("Original DNA Sequence:", dna_seq)
    print("Reverse Complement:", reverse_complement)

    # Transcribe DNA to RNA
    rna_seq = dna_seq.transcribe()
    print("Transcribed mRNA:", rna_seq)

    # Translate RNA to Protein
    protein_seq = rna_seq.translate()
    print("Protein Sequence:", protein_seq)

fasta_file = "/Users/giuse/pythonProject/Mycodes/samples/sequence1.txt"
DNA_translator(fasta_file)
# Execution time
end_time = time.time()
print("Execution Time:", round((end_time - start_time) * 1000, 2), "ms")


from Bio.Seq import Seq
from Bio import SeqIO


def DNA_translate(sequence) -> str:
    """ Load sequences from a list of FASTA files."""
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
    if len(sequence) == 0:
        raise ValueError(f"No sequences found in {fasta_file}.")

    """Translate the nucleotide sequence in amino acids sequence using codon table"""
    if len(sequence) % 3 == 0 :
        protein_sequence = Seq.translate(Seq(sequence), "Standard")
    elif len(sequence) - 1 % 3 == 0:
        protein_sequence = Seq.translate(Seq(sequence)[:-1], "Standard")
    else:
        protein_sequence = Seq.translate(Seq(sequence)[:-2], "Standard")
    return protein_sequence


fasta_file = "/Users/giuse/pythonProject/Mycodes/samples/sequence1.txt" # add sequence file path
sequence = DNA_translate(fasta_file)
print(f"The sequence is {sequence}")
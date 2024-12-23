from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
import subprocess
import os

class FastaFileParser:
    def __init__(self, file_paths: list[str]):
        self.file_paths = file_paths
        self.sequences = self.parse_files()

    def parse_files(self) -> list[SeqRecord]:
        all_sequences = []
        for file_path in self.file_paths:
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"File not found: {file_path}")
            try:
                for record in SeqIO.parse(file_path, "fasta"):
                    all_sequences.append(record)
            except Exception as e:
                raise ValueError(f"Error parsing FASTA file {file_path}: {e}")
        return all_sequences

class MultipleSequenceAligner:
    def __init__(self, sequences: list[SeqRecord], tool: str = "clustalw"):
        self.sequences = sequences
        self.combined_sequences_path = "Combined_Sequences.fasta"
        self.alignment_path = "Combined_Sequences.aln" # Path for the alignment file
        self.tool = tool
        self.format = "clustal" if self.tool == "clustalw" else "fasta"
        SeqIO.write(self.sequences, self.combined_sequences_path, "fasta")

    def perform_msa(self) -> None:
        cmd = self.run_command()
        try:
            subprocess.run(cmd)
        except subprocess.CalledProcessError as e:
            print(f"Error during MSA execution:\n{e.stderr}")
            raise # Re-raise the exception after printing the error

    def read_alignment(self) -> AlignIO.MultipleSeqAlignment | None:
        try:
            return AlignIO.read(self.alignment_path, self.format)
        except FileNotFoundError:
            print("Alignment file not found. Did you run perform_msa()?")
            return None

    def __str__(self):
        alignment = self.read_alignment()
        if alignment:
            return str(alignment)
        return ""

    def run_command(self) -> list[str]:
        match self.tool:
            case "clustalw":
                cmd = ["conda", "run", "clustalw", "-align", f"-infile={self.combined_sequences_path}"]
            case "muscle":
                cmd = ["conda", "run", "muscle", "-align", self.combined_sequences_path, "-output", self.alignment_path]
            case _:
                raise ValueError("Unsupported tool. Use 'clustalw' or 'muscle'.")
        return cmd

if __name__ == "__main__":
    fasta_files = ["/Users/giuse/pythonProject/Mycodes/samples/sequence1.txt",
                   "/Users/giuse/pythonProject/Mycodes/samples/sequence2.txt",
                   "/Users/giuse/pythonProject/Mycodes/samples/sequence3.txt"]

    try:
        parser = FastaFileParser(fasta_files)
        aligner = MultipleSequenceAligner(parser.sequences, tool="clustalw") # Choose tool here
        aligner.perform_msa()
        print(aligner)

        parser = FastaFileParser(fasta_files)
        aligner = MultipleSequenceAligner(parser.sequences, tool="muscle") # Choose tool here
        aligner.perform_msa()
        print(aligner)
    except Exception as e:
        print(f"An error occurred: {e}")
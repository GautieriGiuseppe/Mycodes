from Bio import SeqIO
from Bio import AlignIO
import subprocess

class MultipleSequenceAnalyzer():
    def __init__(self, fasta_file):
        """ Load sequences from a list of FASTA files."""
        self.sequences = []
        for fasta_file in fasta_files:
            for record in SeqIO.parse(fasta_file, "fasta"):
                self.sequences.append(record)
        if len(self.sequences) == 0:
            raise ValueError(f"No sequences found in {self}.")


    def save_combined_fasta(self, output_file="combined_sequences.fasta") -> str:
        """ Save multiple sequences into a single FASTA file for MSA."""
        SeqIO.write(self.sequences, output_file, "fasta")
        return output_file

    def perform_msa(self, tool:str):
        """ Perform MSA using ClustalW or Muscle."""
        command = tool
        cmd = run_command(command)
        subprocess.run(cmd)


    def read_alignment(self, format="clustal"):
        """Read and print the alignment from the alignment file."""
        alignment = AlignIO.read("combined_sequences.aln", format)
        print(alignment)


def run_command(command:str):
    match command:
        case  "clustalw":
            cmd = ["conda", "run", "clustalw", "-align", f"-infile={combined_fasta}"]
        case "muscle":
            cmd = ["conda", "run", "muscle", "-align", f"{combined_fasta}", "-output", "combined_sequences.aln"]
        case _:
            raise ValueError("Unsupported tool. Use 'clustalw' or 'muscle'.")
    return cmd

if __name__ == "__main__":
    # List of FASTA files on local pc
    fasta_files = ["/Users/giuse/pythonProject/Mycodes/samples/sequence1.txt",
                   "/Users/giuse/pythonProject/Mycodes/samples/sequence2.txt",
                   "/Users/giuse/pythonProject/Mycodes/samples/sequence3.txt"] # add path of files here

    # Load and combine sequences into a single file
    sequences = MultipleSequenceAnalyzer(fasta_files)
    combined_fasta = sequences.save_combined_fasta()

    # Perform MSA using selected tool
    alignment = MultipleSequenceAnalyzer(combined_fasta).perform_msa(tool="clustalw")  # or change to 'muscle'
    # Read and print the alignment
    tool = "clustalw"
    MultipleSequenceAnalyzer(alignment).read_alignment(format="clustal" if tool == "clustalw" else "fasta")


import os
import time
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
from scipy.cluster.hierarchy import linkage

class FastaFileParses:
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

# Compute pairwise alignments and generate a distance matrix
def compute_distance_matrix(sequences : list[SeqRecord]) -> np.ndarray:
    num_seqs = len(sequences)
    distance_matrix = np.zeros((num_seqs, num_seqs))

    # Create aligner object
    aligner = PairwiseAligner()
    aligner.mode = 'global' # Needleman-Wunsch
    aligner.match_score, aligner.mismatch_score, aligner.open_gap_score, aligner.extend_gap_score = 1, -1, -2, -1

    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            seq1, seq2 = str(sequences[i].seq), str(sequences[j].seq)
            score = aligner.score(seq1, seq2) # compute alignment score
            max_score = min(len(seq1), len(seq2))
            distance = 1.0 - (score / max_score) # similarity -> distance
            distance_matrix[i, j] = distance_matrix[j, i] = distance
    return distance_matrix

# Generate Guide tree
def generate_guide_tree(distance_matrix : np.ndarray, sequence_ids) -> np.ndarray:
    linkage_matrix = linkage(distance_matrix, method="average")

    return linkage_matrix

# Progressive Alignment
def progressive_alignment(sequences, linkage_matrix):

    # Convert sequences into a Pandas DataFrame
    df = pd.DataFrame({
        "ID": [seq.id for seq in sequences],
        "Sequence": [str(seq.seq) for seq in sequences]
    })

    num_seqs = len(sequences)

    # Create a list to track current alignments
    aligned_sequences = df["Sequence"].values.tolist()
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score, aligner.mismatch_score, aligner.open_gap_score, aligner.extend_gap_score = 1, -1, -2, -1


    # Follow guide tree order for alignment
    for node in linkage_matrix:
        seq1_index, seq2_index = int(node[0]), int(node[1])

        # Extract sequences/clusters to align
        seq1 = aligned_sequences[seq1_index]
        seq2 = aligned_sequences[seq2_index]

        # Perform pairwise alignment
        alignment = aligner.align(seq1, seq2)[0] # Get best alignment
        aligned_seq1, aligned_seq2 = alignment.aligned # extract alignmed sequences

        # Store aligned sequences
        aligned_sequences.append(str(alignment))

    # Update DataFrame
    df["Aligned_Sequence"] = aligned_sequences[:num_seqs]

    return df


# Usage
start_time = time.time()
if __name__ == '__main__':

    fasta_files = ["/Users/giuse/pythonProject/Mycodes/samples/sequence1.txt",
                   "/Users/giuse/pythonProject/Mycodes/samples/sequence2.txt",
                   "/Users/giuse/pythonProject/Mycodes/samples/sequence3.txt"]
    parser = FastaFileParses(fasta_files)
    sequences = parser.sequences

    if len(sequences) < 2:
        raise ValueError("At least two sequences required.")

    sequence_ids = [record.id for record in sequences]

    # Compute distance matrix
    distance_matrix = compute_distance_matrix(sequences)

    # Generate Guide Tree
    guide_tree = generate_guide_tree(distance_matrix, sequence_ids)

    # Perform progressive alignment
    msa_df = progressive_alignment(sequences, guide_tree)

    print(msa_df[["ID", "Aligned_Sequence"]])

    end_time = time.time()
    print(f"The execution time is {end_time - start_time} seconds")
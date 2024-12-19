import random

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1):
    # Create score matrix and traceback matrix
    m, n = len(seq1), len(seq2)
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]
    traceback_matrix = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize the score matrix with gap penalties
    for i in range(1, m + 1):
        score_matrix[i][0] = gap * i
    for j in range(1, n + 1):
        score_matrix[0][j] = gap * j

    # Fill the score matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diag = score_matrix[i - 1][j - 1] + match_score
            up = score_matrix[i - 1][j] + gap
            left = score_matrix[i][j - 1] + gap
            score_matrix[i][j] = max(diag, up, left)

    # Traceback to get alignment
    aligned_seq1, aligned_seq2 = '', ''
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i - 1][j - 1] + (
        match if seq1[i - 1] == seq2[j - 1] else mismatch):
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif i > 0 and score_matrix[i][j] == score_matrix[i - 1][j] + gap:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1

    return aligned_seq1, aligned_seq2


# Example usage
bases = 'ATGC'
sequence = ''.join([bases[random.randint(0,3)]for x in range(1000)])


seq1 = ''.join([bases[random.randint(0,3)]for x in range(1000)])
seq2 = ''.join([bases[random.randint(0,3)]for x in range(1000)])
seq3 = ''.join([bases[random.randint(0, 3)]for x in range(1000)])


aligned_seq1, aligned_seq2 = needleman_wunsch(seq1, seq2)
aligned_seq2, aligned_seq3 = needleman_wunsch(seq2, seq3)
print("Aligned Sequence 1:", aligned_seq1)
print("Aligned Sequence 2:", aligned_seq2)
print("Aligned Sequence 3:", aligned_seq3)

# MSA basic approach
def progressive_alignment(sequences):
    # Start with the first sequence
    aligned_seq = sequences[0]

    for i in range(1, len(sequences)):
        aligned_seq, _ = needleman_wunsch(aligned_seq, sequences[i])

    return aligned_seq


# Example usage for multiple sequences
sequences = [aligned_seq1, aligned_seq2, "AGTACAGTTC"]
aligned = progressive_alignment(sequences)
print("alignment result:  ", aligned)
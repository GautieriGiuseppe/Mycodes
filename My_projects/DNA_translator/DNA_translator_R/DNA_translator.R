# DNA Sequence Translator
library(Biostrings)

start_time <- Sys.time()
# Load FASTA
fasta_file <- "/Users/giuse/pythonProject/Mycodes/samples/sequence1.txt"
dna_sequence <- readDNAStringSet(fasta_file)

# Print Original DNA Sequence
cat("Original DNA sequence:\n")
print(dna_sequence)

# Get Reverse Complement
reverse_complement <- reverseComplement(dna_sequence)
cat("\nReverse Complement:\n")
print(reverse_complement)

# Transcribe to mRNA (DNA -> RNA)
rna_sequence <- RNAStringSet(dna_sequence)
cat("\nTranscribed mRNA:\n")
print(rna_sequence)

# Translate (RNA -> Protein)
protein_sequence <- translate(dna_sequence)
cat("\nProtein Sequence:\n")
print(protein_sequence)

end_time <- Sys.time()
exec_time <- end_time - start_time
print(exec_time * 1000) # milliseconds
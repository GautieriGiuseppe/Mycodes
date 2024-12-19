import re
from Bio import SeqIO


class Sequence():
    def __init__(self, fasta_file):
        """ Load sequences from a list of FASTA files."""
        for record in SeqIO.parse(fasta_file, "fasta"):
            self.sequence = str(record.seq)
            self.id = record.id
            self.description = record.description
            self.length = len(self.sequence)
            self.results = []


    def calculate_gc_content(self):
        """Calculate the GC content of a DNA sequence."""
        g = self.sequence.count('G')
        c = self.sequence.count('C')
        return (g + c) / self.length * 100


    def calculate_nucleotide_counts_DNA(self):
        """Calculate the count of each nucleotide in a DNA sequence."""
        return {
            'A': self.sequence.count('A'),
            'T': self.sequence.count('T'),
            'G': self.sequence.count('G'),
            'C': self.sequence.count('C')
        }


    def calculate_nucleotide_counts_RNA(self):
        """Calculate the count of each nucleotide in a DNA sequence."""
        return {
            'A': self.sequence.count('A'),
            'U': self.sequence.count('U'),
            'G': self.sequence.count('G'),
            'C': self.sequence.count('C')
        }


    def find_promoter_sequences(self, promoter_patterns):
        """Identify promoter sequences in a DNA sequence based on given patterns."""
        promoters = []
        for pattern in promoter_patterns:
            for match in re.finditer(pattern, self.sequence):
                promoters.append((match.start(), match.end(), match.group()))
        return promoters


    def find_transcription_start_site(self, promoters):
        """Identify the transcription start site using promoter distance and sequence."""
        start_site_pattern = '[CT][CT]A[ATGC][AT][CT][CT]'
        start_site = []
        distance = 0
        for match in re.finditer(start_site_pattern, self.sequence):
            start_site.append((match.start(), match.end(), match.group()))
        for start_init, end_init, init in start_site:
            for start_prom, end_prom, promoter in promoters:
                distance = abs(start_init - end_prom)
        return start_site, distance


    def analyze_genome(self, promoter_patterns):
        """Analyze a genome from a FASTA file, including GC content, nucleotide counts, and promoter sequences."""
        #for record in SeqIO.parse(self, "fasta"):
        gc_content = self.calculate_gc_content()
        nucleotide_counts = self.calculate_nucleotide_counts_DNA()
        promoters = self.find_promoter_sequences(promoter_patterns)
        (start_site, distance) = self.find_transcription_start_site(promoters)
        self.results.append({
            'id': self.id,
            'description': self.description,
            'length': self.length,
            'gc_content': gc_content,
            'nucleotide_counts': nucleotide_counts,
            'promoters': promoters,
            'transcription_start_site': start_site,
            'distance_promoter': distance,
        })
        return self.__str__()

    def __str__(self):
        """Displays the results"""
        for result in self.results:
            print(f"ID: {result['id']}")
            print(f"Description: {result['description']}")
            print(f"Length: {result['length']}")
            print(f"GC Content: {result['gc_content']:.2f}%")
            print("Nucleotide Counts:")
            for nucleotide, count in result['nucleotide_counts'].items():
                print(f"  {nucleotide}: {count}")
            print("Promoter Sequences:")
            for start, end, promoter in result['promoters']:
                print(f"  Position: {start}-{end}, Sequence: {promoter}")
            print("Transcription start site: ")
            if result['transcription_start_site']:
                for start, end, start_site in result['transcription_start_site']:
                    print(f"  Position: {start}-{end}, Sequence: {start_site}")
                    if result['distance_promoter'] < 25 or result['distance_promoter'] > 30:
                        print(f"  Distance: {result['distance_promoter']}, Abnormal promoter position")
                    else:
                        print(f"  Distance: {result['distance_promoter']}, Physiological positioning")
            else:
                print()
                print("  Transcription start site not found in this sequence")
            print()

fasta_file = "" # add path to fasta file here
sequence = Sequence(fasta_file)
promoter_patterns = [r'TATA[AT]A[AT]', r'CAAT', r'CGTACG']  # Example patterns: TATA box, CAAT box, etc.
analysis_results = sequence.analyze_genome(promoter_patterns) # here we call the main method


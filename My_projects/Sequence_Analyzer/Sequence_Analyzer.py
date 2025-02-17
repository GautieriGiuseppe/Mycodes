import os
import re
from typing import *
from Bio import SeqIO

class ShortSeq:
    def __init__(self, startPos: int, endPos: int, seq: str) -> None:
        self.startPos, self.endPos, self.seq = startPos, endPos, seq

    def fromPattern(pattern: Pattern[str] | str, seq: str) -> list:
        return [ShortSeq(m.start(), m.end(), m.group()) for m in re.finditer(pattern, seq)]

    def __repr__(self) -> str:
        return f"Position: {self.startPos}-{self.endPos}, Sequence: {self.seq}"

class Sequence:
    PROMOTER_PATTERNS = tuple(map(re.compile, (r'TATA[AT]A[AT]', r'CAAT', r'CGTACG')))

    def fromFasta(filePath: str) -> list:
        if not os.path.exists(filePath):
            raise FileNotFoundError(f"File not found: {filePath}")
        try:
            return [Sequence.fromFastaRecord(record) for record in SeqIO.parse(filePath, "fasta")]
        except Exception as e:
            raise ValueError(f"Error parsing FASTA file: {e}")

    def fromFastaRecord(record) -> "Sequence":
        return Sequence(str(record.seq), record.id, record.description)

    def __init__(self, seq: str, id=0, descr="Not provided.", promoterPatterns: tuple[Pattern[str], ...] = PROMOTER_PATTERNS) -> None:
        self.id, self.sequence, self.description = id, seq, descr
        self.gc_content = self.calculate_gc_content()
        self.promoters = self.findAllShortSeqsFromPatterns(promoterPatterns)
        self.nucleotide_counts = self.calculate_nucleotide_counts_DNA()
        self.start_sites, self.dists = self.find_transcription_start_site(self.promoters)

    def calculate_gc_content(self) -> float:
        g = self.sequence.count('G')
        c = self.sequence.count('C')
        return (g + c) / len(self.sequence) * 100

    def calculate_nucleotide_counts_DNA(self) -> dict[str, int]:
        return {
            'A': self.sequence.count('A'),
            'T': self.sequence.count('T'),
            'G': self.sequence.count('G'),
            'C': self.sequence.count('C')
        }

    def findAllShortSeqsFromPatterns(self, patterns: tuple[Pattern[str] | str, ...]) -> list[ShortSeq]:
        return [seq for pattern in patterns for seq in ShortSeq.fromPattern(pattern, self.sequence)]

    def find_transcription_start_site(self, promoters: list[ShortSeq]) -> tuple[list[ShortSeq], list[int]]:
        startSites = self.findAllShortSeqsFromPatterns([r'[CT][CT]A[ATGC][AT][CT][CT]'])
        distsToPromoters = [
            abs(transcrSeq.startPos - promoter.endPos)
            for transcrSeq in startSites
            for promoter in promoters
        ]
        return startSites, distsToPromoters

    def __str__(self):
        seqRepr = (
            f"ID: {self.id}\nDescription: {self.description}\nLength: {len(self.sequence)}\nGC Content: {self.gc_content:.2f}%\n\n"
            f"Nucleotide Counts: {self.nucleotide_counts}\n"
        )
        if self.start_sites:
            seqRepr += "Transcription Start Sites:\n"
            for i, transcrSeq in enumerate(self.start_sites):
                d = self.dists[i] if i < len(self.dists) else "N/A"
                promoter_position = "Abnormal" if d == "N/A" or d < 25 or d > 30 else "Physiological"
                seqRepr += f"{transcrSeq}, Distance: {d}, {promoter_position}\n"
        else:
            seqRepr += "Transcription Start Sites: not found\n"
        return seqRepr

if __name__ == "__main__":
    try:
        sequences = Sequence.fromFasta("")  # Update this to your file path
        for seq in sequences:
            print(seq)
    except Exception as e:
        print(f"Error: {e}")
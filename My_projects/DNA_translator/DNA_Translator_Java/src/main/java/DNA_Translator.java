import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.template.SequenceView;
import org.biojava.nbio.core.sequence.transcription.TranscriptionEngine;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.RNASequence;

import java.io.File;
import java.util.LinkedHashMap;

public class DNA_Translator {
    public static void main(String[] args) {
        long startTime = System.nanoTime();

        try {
            // Load DNA sequence from FASTA file
            File fastaFile = new File("/Users/giuse/pythonProject/Mycodes/samples/sequence1.txt");
            LinkedHashMap<String, DNASequence> sequences = FastaReaderHelper.readFastaDNASequence(fastaFile);

            // Get the sequence
            DNASequence dna = sequences.values().iterator().next();
            System.out.println("Loaded DNA Sequence: " + dna.getSequenceAsString());

            // Get Reverse Complement
            SequenceView reverseComplementView = dna.getReverseComplement();
            DNASequence reverseComplement = new DNASequence(reverseComplementView.getSequenceAsString());
            System.out.println("Original DNA Sequence: " + dna.getSequenceAsString());
            System.out.println("Reverse Complement: " + reverseComplement.getSequenceAsString());


            // Transcribe DNA to RNA
            TranscriptionEngine engine = TranscriptionEngine.getDefault();
            RNASequence rna = dna.getRNASequence(engine);
            System.out.println("Transcribed mRNA: " + rna.getSequenceAsString());

            // Translate RNA to Protein
            ProteinSequence protein = rna.getProteinSequence(engine);
            System.out.println("Protein Sequence: " + protein.getSequenceAsString());

        } catch (Exception e) {
            e.printStackTrace();
        }
        long endTime = System.nanoTime();
        long duration = endTime - startTime;
        System.out.println("Execution Time " + duration / 1_000_000 + "ms");
    }
}
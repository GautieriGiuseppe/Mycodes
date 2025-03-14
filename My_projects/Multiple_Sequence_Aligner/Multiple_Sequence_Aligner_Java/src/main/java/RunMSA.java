import java.io.*;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.*;

public class RunMSA {

    public static void main(String[] args){
        List<String> fastaFiles = Arrays.asList("/Users/giuse/pythonProject/Mycodes/samples/sequence1.txt",
                                                "/Users/giuse/pythonProject/Mycodes/samples/sequence2.txt",
                                                "/Users/giuse/pythonProject/Mycodes/samples/sequence3.txt");

        // Define output files for ClustalW and MUSCLE
        String combinedFasta = "combined.fasta";
        String clustalwOutput = "aligned_clustalw.fasta";
        String muscleOutput = "aligned_muscle.fasta";

        // Start measure time
        long startCombine = System.nanoTime();
        combineFastaFiles(fastaFiles, combinedFasta);
        long endCombine = System.nanoTime();
        System.out.println("FASTA Merging Time: " + ((endCombine - startCombine) / 1e9) + " seconds");

        // Create ExecutorService to run ClustalW and MUSCLE in parallel
        long startTime = System.nanoTime();
        ExecutorService executor = Executors.newFixedThreadPool(2);

        Future<Double> clustalwTime = executor.submit(() -> runClustalW(combinedFasta, clustalwOutput));
        Future<Double> muscleTime = executor.submit(() -> runMuscle(combinedFasta, muscleOutput));

        try {
            System.out.println("\nExecution Times: ");
            System.out.println("ClustalW Time: " + clustalwTime.get() + " seconds");
            System.out.println("Muscle Time: " + muscleTime.get() + " seconds");

            // Read & Print alignment results
            System.out.println("\n*********** ClustalW Alignment Result ***********");
            printFileContents(clustalwOutput);
            System.out.println("\n*********** Muscle Alignment Result ***********");
            printFileContents(muscleOutput);

        } catch (Exception e) {
            e.printStackTrace();
        }

        executor.shutdown();
        long endTime = System.nanoTime();
        double totalExecutionTime = (endTime - startTime) / 1e9;
        System.out.println("Total Execution Time: " + totalExecutionTime + " seconds");
    }

    // Merge FASTA files in a single combined.fasta
    public static void combineFastaFiles(List<String> fastaFiles, String outputFile) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
            for (String fastaFile : fastaFiles) {
                try (BufferedReader reader = new BufferedReader(new FileReader(fastaFile))) {
                    String line;
                    while ((line = reader.readLine()) != null) {
                        writer.write(line);
                        writer.newLine(); // Maintain FASTA format
                    }
                }
            }
        } catch (IOException e){
            e.printStackTrace();
        }
    }

    // Run ClustalW with multiple FASTA files
    public static double runClustalW(String inputFile, String outputFile){
        long startTime = System.nanoTime();
        try {
            ProcessBuilder pb = new ProcessBuilder(
                    "conda", "run", "clustalw", "-align", "-INFILE=" + inputFile, "-OUTFILE=" + outputFile);
            Process process = pb.start();
            captureProcessOutput(process);
            process.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
        long endTime = System.nanoTime();
        return (endTime - startTime) / 1e9;
    }

    // Run MUSCLE with multiple FASTA files
    public static double runMuscle(String inputFile, String outputFile) {
        long startTime = System.nanoTime();
        try {
            ProcessBuilder pb = new ProcessBuilder("conda", "run", "muscle", "-align", inputFile, "-output", outputFile);
            Process process = pb.start();
            captureProcessOutput(process);
            process.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
        long endTime = System.nanoTime();
        return (endTime - startTime) / 1e9;
    }

    // Handle output and errors
    private static void captureProcessOutput(Process process){
        try (BufferedReader stdInput = new BufferedReader(new InputStreamReader(process.getInputStream()));
            BufferedReader stdError = new BufferedReader(new InputStreamReader(process.getErrorStream()))){

            String line;
            while((line = stdInput.readLine()) != null){
                System.out.println(line);
            }
            while((line = stdError.readLine()) != null){
                System.err.println(line);
            }
        } catch (IOException e){
            e.printStackTrace();
        }
    }

    // Read and Print file contents
    private static void printFileContents(String filename){
        try (BufferedReader reader = new BufferedReader(new FileReader(filename))){
            String line;
            while ((line = reader.readLine()) != null){
                System.out.println(line);
            }
        } catch (IOException e){
            System.err.println("[ERROR] Cannot read file: " + filename);
        }
    }
}


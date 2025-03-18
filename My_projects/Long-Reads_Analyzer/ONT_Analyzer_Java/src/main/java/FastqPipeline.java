import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.*;

public class FastqPipeline {
    private final String basePath;
    private final String resultsPath;
    private final String inputPath;
    private final String referenceGenome;
    private final String summaryFile;
    private final String condaEnv;
    private final int numWorkers;
    private final List<Integer> barcodes = Arrays.asList(1, 2, 3, 4, 5, 6); // Barcodes to process

    public FastqPipeline(String basePath, String referenceGenome, String summaryFile, int numWorkers){
        this.basePath = basePath;
        this.resultsPath = basePath + "/results";
        this.inputPath = basePath + "/data_separated/";
        this.referenceGenome = referenceGenome;
        this.summaryFile = summaryFile;
        this.condaEnv = "long_reads";
        this.numWorkers = (numWorkers > 0) ? numWorkers : Runtime.getRuntime().availableProcessors();

        // Ensure results directory exists
        new File(resultsPath).mkdirs();
    }

    private int runCommand(List<String> command){
        try {
            ProcessBuilder pb = new ProcessBuilder(command);
            pb.redirectErrorStream(true);
            Process process = pb.start();

            try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))){
                String line;
                while ((line = reader.readLine()) != null){
                    System.out.println(line);
                }
            }
            return process.waitFor();
        } catch (IOException | InterruptedException e){
            e.printStackTrace();
            return 1;
        }
    }

    public void runQualityCheck(){
        System.out.println("Running PycoQC...");
        List<String> cmd = Arrays.asList("conda", "run", "-n", condaEnv, "pycoQC", "--summary_file", summaryFile,
                                        "--html_outfile", resultsPath + "/Run_quality.html");
        runCommand(cmd);
    }

    public void processBarcode(int barcode, String tool){
        List<String> cmd = new ArrayList<>();
        String inputFile, outputFile;

        switch (tool){
            case "filtlong":
                inputFile = inputPath + "barcode_" + barcode + ".fastq";
                outputFile = resultsPath + "/barcode_" + barcode + ".filtered.fastq";
                cmd = Arrays.asList("conda", "run", "filtlong", "--min_length", "1000", "--keep_percent", "90",
                                    "--trim", "-a", referenceGenome, inputFile);
                break;

            case "fastqc":
                inputFile = resultsPath + "/barcode_" + barcode + ".filtered.fastq";
                cmd = Arrays.asList("conda", "run", "-n", condaEnv, "fastqc", inputFile, "--outdir", resultsPath + "/fastqc");
                break;
            case "mapping":
                inputFile = resultsPath + "/barcode_" + barcode + ".filtered.fastq";
                String samOutput = resultsPath + "/sample_" + barcode + ".sam";
                cmd = Arrays.asList("minimap2", "-ax", "map-ont", "-t", "4", referenceGenome, inputFile, ">", samOutput);
                break;
            case "samtools":
                String inputBam = resultsPath + "/sample_" + barcode + ".bam";
                cmd = Arrays.asList("samtools", "view", "-bS", "-@", "4", resultsPath + "/sample_" + barcode + ".sam",
                                    "|", "samtools", "sort", "-o", inputBam);
                break;
            case "picard":
                inputBam = resultsPath + "/sample_" + barcode + ".bam";
                outputFile = resultsPath + "/sample_" + barcode + ".markdup.bam";
                cmd = Arrays.asList("conda", "run", "-n", condaEnv, "picard", "_Xmx8G", "MarkDuplicates",
                                    "I=" + inputBam, "O=" + outputFile, "REMOVE_DUPLICATES=true", "CREATE_INDEX=true");
                break;
            case "variant_calling":
                inputBam = resultsPath + "/sample_" + barcode + "markdup.bam";
                outputFile = resultsPath + "/sample_" + barcode + ".vcf";
                cmd = Arrays.asList("conda", "run", "-n", condaEnv, "sniffles", "--input", inputBam, "--vcf", outputFile);
                break;
            default:
                throw new IllegalArgumentException("Unsupported tool: " + tool);
        }

        System.out.println("Running " + tool + " on barcode " + barcode + "...");
        runCommand(cmd);
    }

    public void runParallel(String tool){
        System.out.println("Starting parallel processing for: " + tool);
        ExecutorService executor = Executors.newFixedThreadPool(numWorkers);
        List<Future<?>> futures = new ArrayList<>();

        for (int barcode : barcodes){
            int finalBarcode = barcode;
            futures.add(executor.submit(() -> processBarcode(finalBarcode, tool)));

        }

        for (Future<?> future : futures){
            try {
                future.get(); // wait all tasks to complete
            } catch (InterruptedException | ExecutionException e){
                e.printStackTrace();
            }
        }

        executor.shutdown();
    }

    public void runPipeline(){
        System.out.println("Starting ONT sequencing pipeline...");

        // Step 1: Quality check of the run
        runQualityCheck();
        System.out.println("Run Quality Check Completed!");

        // Step 2: Parallel Processing
        runParallel("filtlong");
        runParallel("fastqc");
        runParallel("mapping");
        runParallel("samtools");
        runParallel("picard");
        runParallel("variant_calling");

        System.out.println("Pipeline completed successfully!");
    }

    public static void main(String[] args) {
        String BASE_PATH = "/Users/giuse/pythonProject/Mycodes/MyProjects";
        String REFERENCE_GENOME_PATH = BASE_PATH + "/Escherichia_coli_reference";
        String SUMMARY_FILE_PATH = BASE_PATH + "/ONT_simulated_summary.txt";

        // Check files
        for (String path : Arrays.asList(REFERENCE_GENOME_PATH, SUMMARY_FILE_PATH)){
            if (!Files.exists(Paths.get(path))){
                throw new RuntimeException("File not found: " + path);
            }
        }

        // Run
        FastqPipeline pipeline = new FastqPipeline(BASE_PATH, REFERENCE_GENOME_PATH, SUMMARY_FILE_PATH, 4);
        pipeline.runPipeline();
    }
}

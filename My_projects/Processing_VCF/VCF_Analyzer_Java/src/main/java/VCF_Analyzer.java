import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.csv.*;



import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;

public class VCF_Analyzer {
    private List<String> inputVcfs;
    private String mergedVcf = "merged_output.vcf";
    private String clinvarFile;
    private String geneDbFile;
    private List<Map<String, String>> variantData = new ArrayList<>();

    public VCF_Analyzer(List<String> inputVcfs, String clinvarFile, String geneDbFile){
        this.inputVcfs = inputVcfs;
        this.clinvarFile = clinvarFile;
        this.geneDbFile = geneDbFile;
    }

    // Merge Multiple VCFs into one
    public void mergeVcfs() throws IOException{
        if (inputVcfs.size() == 1){
            System.out.println("Only one found, skipping merging.");
            mergedVcf = inputVcfs.get(0);
            return;
        }

        System.out.println("Merging VCF files...");
        try (BufferedWriter writer = Files.newBufferedWriter(Paths.get(mergedVcf))){
            for (String vcf : inputVcfs){
                try (BufferedReader reader = Files.newBufferedReader(Paths.get(vcf))){
                    String line;
                    while ((line = reader.readLine()) != null){
                        if (!line.startsWith(("#")) || vcf.equals(inputVcfs.get(0))){
                            writer.write(line);
                            writer.newLine();
                        }
                    }
                }
            }
        }
        System.out.println("Merged VCF saved to " + mergedVcf);


    }

    // Convert VCF to Java List
    public void vcfToListParallel(int numThreads){
        System.out.println("Converting VCF to List in parallel...");
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);
        List<Future<List<Map<String, String>>>> futures = new ArrayList<>();

        try (VCFFileReader vcfReader = new VCFFileReader(new File(mergedVcf), false)){
            VCFHeader header = vcfReader.getFileHeader();
            List<VariantContext> variants = new ArrayList<>();
            vcfReader.iterator().forEachRemaining(variants::add);

            int chunkSize = 50000;
            for (int i = 0; i < variants.size(); i += chunkSize){
                List<VariantContext> chunk = variants.subList(i, Math.min(i + chunkSize, variants.size()));
                futures.add(executor.submit(() -> processVariantChunk(chunk)));
            }

            for (Future<List<Map<String, String>>> future : futures){
                variantData.addAll(future.get());
            }
        } catch (Exception e) {
            e.printStackTrace();
        }finally {
            executor.shutdown();
        }

        System.out.println("Converted " + variantData.size() + " variants into List.");

    }

    private List<Map<String, String>> processVariantChunk(List<VariantContext> chunk){
        List<Map<String, String>> processedData = new ArrayList<>();
        for (VariantContext variant : chunk){
            Map<String, String> row = new HashMap<>();
            row.put("CHROM", variant.getContig());
            row.put("POS", String.valueOf(variant.getStart()));
            row.put("ID", variant.getID());
            row.put("REF", variant.getReference().getBaseString());
            row.put("ALT", variant.getAlternateAlleles().toString());
            row.put("QUAL", String.valueOf(variant.getPhredScaledQual()));
            row.put("FILTER", variant.getFilters().toString());
            processedData.add(row);
        }
        return processedData;
    }

    // Filter Variants
    public void filterVariants(double minQual){
        System.out.println("Filtering variants...");
        int beforeCount = variantData.size();
        variantData = variantData.stream()
                .filter(v -> Double.parseDouble(v.get("QUAL")) >= minQual)
                .collect(Collectors.toList());
        System.out.println("Remaining variants after filtering: " + variantData.size() + " (removed " + (beforeCount - variantData.size()) + ")");
    }

    // Count SNPs and Indels
    public void countVariantTypes(){
        System.out.println("Counting SNPs and Indels...");
        Map<String, Long> counts = variantData.stream()
                .collect(Collectors.groupingBy(v -> v.get("REF").length() == 1 && v.get("ALT").length() == 1 ? "SNP" : "INDEL", Collectors.counting()));

        counts.forEach((key, value) -> System.out.println(key + ": " + value));
    }

    // Find ClinVar Matches
    public void findClinVarVariants() throws IOException{
        System.out.println("Searching for ClinVar variants...");
        Set<String> clinvarPositions = new HashSet<>();

        try (Reader reader = Files.newBufferedReader(Paths.get(clinvarFile));
            CSVParser csvParser = new CSVParser(reader, CSVFormat.DEFAULT.withDelimiter('\t').withFirstRecordAsHeader())){
            for (CSVRecord record : csvParser){
                clinvarPositions.add(record.get("CHROM") + ":" + record.get("POS"));

            }
        }

        long matches = variantData.stream()
                .filter(v -> clinvarPositions.contains(v.get("CHROM") + ":" + v.get("POS")))
                .count();
        System.out.println("Found " + matches + " ClinVar matches.");
    }

    // Annotate Variants
    public void annotateVariants() throws IOException{
        System.out.println("Annotating variants...");
        Map<String, String> geneDict = new HashMap<>();

        try (Reader reader = Files.newBufferedReader(Paths.get(geneDbFile));
             CSVParser csvParser = new CSVParser(reader, CSVFormat.TDF.withFirstRecordAsHeader())) {
            for (CSVRecord record : csvParser){
                geneDict.put(record.get("CHROM") + ":" + record.get("POS"), record.get("GENE"));
            }
        }

        variantData.forEach(v -> v.put("GENE", geneDict.getOrDefault(v.get("CHROM") + ":" + v.get("POS"), "Unknown")));
        System.out.println("Annotation completed.");
    }

    // Save to CSV
    public void saveDataFrame() throws IOException{
        System.out.println("Saving to CSV...");
        try (BufferedWriter writer = Files.newBufferedWriter(Paths.get("processed_variants.csv"));
             CSVPrinter csvPrinter = new CSVPrinter(writer, CSVFormat.DEFAULT.withHeader("CHROM", "POS", "ID", "REF","ALT", "QUAL", "FILTER", "GENE"))){
            for (Map<String, String> row : variantData) {
                csvPrinter.printRecord(row.values());
            }
        }
        System.out.println("Saved processed variants to CSV.");
    }

    public static void main(String[] args) throws IOException{
        VCF_Analyzer analyzer = new VCF_Analyzer(Arrays.asList("sample1.vcf", "sample2.vcf"), "clinvar.tsv", "gene_database.tsv");
        analyzer.mergeVcfs();
        analyzer.vcfToListParallel(Runtime.getRuntime().availableProcessors());
        analyzer.filterVariants(30.0);
        analyzer.countVariantTypes();
        analyzer.findClinVarVariants();
        analyzer.annotateVariants();
        analyzer.saveDataFrame();
        System.out.println("VCF processing completed.");
    }
}

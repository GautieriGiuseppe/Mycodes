from Bio import SeqIO
import subprocess


class Long_Read(str):
    def __init__(self, file):
        """Load and parse the fastq file"""
        for record in SeqIO.parse(file, "fastq"):
            self.name = record.name
            self.sequence = record.seq
            self.id = record.id

    def run_quality_check(self, summary_file):
        """Perform run quality check using PycoQC"""
        act = ["conda", "activate", "long_reads"]
        subprocess.run(act, shell=True)
        cmd = ["pycoQC", f"{summary_file}", "-o", "pycoQC_out.html"]
        subprocess.run(cmd)
        return "Run Quality Check completed"

    def demultiplex(self, tool = "guppy_barcoder"):
        """Demultiplex to separate the samples using guppy or porechop"""
        if tool == "guppy_barcoder":
            cmd = ["bin/guppy_barcoder", "-i", "bin/guppy_out/pass", "-s",
               "bin/guppy_out_dem_trim", "-trim_barcodes"]
        elif tool == "porechop":
            cmd = ["porechop", "-i", "fastq_pass", "-b", "porechop_demultiplexed"]
        else:
            raise ValueError
        subprocess.run(cmd, shell=True)
        return "Reads separated correctly"

    def merge(self):
        """Merge all the fastq files per sample using cat"""
        pass

    def filter_trim_reads(self, rCRS_fasta): # set the parameters of the filtering
        """Filter the reads and trimming using FiltLong."""
        cmd = ["conda", "run", "filtlong", "--min_length", "1000", "--keep_percent 90",
            "--trim", "-1", f"{rCRS_fasta}", f"{self}", ">", "barcode*.trim.fastq"]
        subprocess.run(cmd, shell=True)
        return "Filtering completed"

    def reads_quality_check(self):
        """ Perform quality check using FastQC."""
        cmd = ["conda", "run", "fastqc", f"{self}"] # pass the fastq file
        subprocess.run(cmd)
        return "Quality check completed and available in HTML file"

    def genome_assembly(self, tool="Minimap2"):
        """Perform genome assembly using minimap"""
        if tool == "Flye":
            cmd = ["conda", "run", "flye", "-nano-hq", f"{self}", "-genome-size",
                   "-out-dir", "./flye_assembly/hq", "-threads 7"] # customize the expected genome size
        elif tool == "Minimap2":
            cmd = ["conda", "run", "bin/minimap2", "bin/out/all_reads.fastq",
                   ">", "bin/alignment.sam", "-ax", "map-ont"] # pass fastq with all the reads
        else:
            raise ValueError
        subprocess.run(cmd, shell=True)
        cmd = ["conda", "run", "samtools", "fixmate"]
        subprocess.run(cmd)
        return "Genome Assembly completed"

    def remove_duplicates(self):
        """Remove duplicates using picard"""
        cmd = ["picard", "MarkDuplicates", "I=barcode*.bam", "O=barcode*.markdup.bam",
               "METRICS_FILE=rCRS.metric", "REMOVE_DUPLICATES=true", "CREATE_INDEX=true"]
        subprocess.run(cmd, shell=True)
        return "Duplicates removed successfully"

    def variant_calling(self, tool="cuteSV"):
        """Variant calling using CuteSV and Sniffles"""
        if tool == "cuteSV":
            cmd = ["conda", "run", "cuteSV", "input.markbam", "rCRS.fasta", "output.outcute",
                   "folder in which save output"]
        elif tool == "Sniffles":
            cmd = []
        else:
            raise ValueError
        subprocess.run(cmd, shell=True)
        return "Variant Calling Performed"

    def complete_analysis(self, summary_file, rCRS):
        """Method to perform the complete pipeline"""
        run_quality = self.run_quality_check(summary_file)
        demultiplex = self.demultiplex()
        merge = self.merge() # to complete
        reads_filter = self.filter_trim_reads(rCRS)
        reads_quality_check = self.reads_quality_check()
        genome = self.genome_assembly()
        remove_duplicates = self.remove_duplicates()
        variant_calling = self.variant_calling()

    def __str__(self):
        pass



fastq_file = Long_Read("/Users/giuse/pythonProject/Mycodes/samples/1_control_ITS2_2019_minq7.fastq")
fastq_file.complete_analysis()
fastq_file.run_quality_check(summary_file="")

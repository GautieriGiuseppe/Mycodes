import subprocess
import os
import Run_Quality_Analyzer

""" Example fastq composition per sequence, many in the same file
@SRR077391.1 HWUSI-EAS667_105020215:3:1:1339:1030/2
TCCCTTTCTACTTAAATCTTGGTTGAGTTTAGGTGTCATTGTTTTGTATACCTACTGAATCCATATTTGTCTGCATCAGAGCACATTATGACATCTGGGT
+
G3BGBGE@===?4BB<>59ABCA=4-=B?B<EB:D8+,==BBA87+<+@6BD8-C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""


def run_quality_check(summary_file:str) -> subprocess.CompletedProcess:
    """Perform run quality check using PycoQC"""
    Run_Quality_Analyzer.Graphs(summary_file)
    cmd = ["conda", "run", "-n", "long_reads", "pycoQC", "--summary_file", summary_file]
    return subprocess.run(cmd, text=True, check=True)

class Fastq:
    def __init__(self, file: str, reference_genome: str = None) -> None:
        """Load and parse the FASTQ compressed file"""
        self.file = file
        self.reference_genome = reference_genome

    def demultiplex(self, tool="guppy_barcoder") -> subprocess.CompletedProcess:
        """Demultiplex samples using Guppy or Porechop"""
        cmd = self.run_command(tool)
        return subprocess.run(cmd, text=True, check=True)

    def filter_trim_reads(self) -> subprocess.CompletedProcess:
        """Filter and trim reads using FiltLong"""
        for sample in range(1, 11):
            cmd = self.run_command("filtlong", sample)
            subprocess.run(cmd, check=True, text=True)
        return subprocess.CompletedProcess(args=cmd, returncode=0)

    def reads_quality_check(self) -> subprocess.CompletedProcess:
        """Perform quality check using FastQC"""
        for sample in range(1, 11):
            cmd = self.run_command("fastqc", sample)
            subprocess.run(cmd, check=True)
        return subprocess.CompletedProcess(args=cmd, returncode=0)

    def genome_assembly(self, tool="Minimap2") -> subprocess.CompletedProcess:
        """Perform genome assembly"""
        for sample in range(1, 11):
            cmd = self.run_command(tool, sample)
            subprocess.run(cmd, check=True)
        return subprocess.CompletedProcess(args=cmd, returncode=0)

    def sort_index(self) -> subprocess.CompletedProcess:
        """Sort and index using samtools"""
        for sample in range(1, 11):
            cmd = self.run_command("samtools_view", sample)
            subprocess.run(cmd, check=True)
            cmd = self.run_command("samtools_index", sample)
            subprocess.run(cmd, check=True)
        return subprocess.CompletedProcess(args=cmd, returncode=0)

    def remove_duplicates(self) -> subprocess.CompletedProcess:
        """Remove duplicates using Picard"""
        for sample in range(1, 11):
            cmd = self.run_command("picard", sample)
            subprocess.run(cmd, check=True)
        return subprocess.CompletedProcess(args=cmd, returncode=0)

    def variant_calling(self, tool="cuteSV") -> subprocess.CompletedProcess:
        """Perform variant calling using CuteSV or Sniffles"""
        for sample in range(1, 11):
            cmd = self.run_command(tool, sample)
            subprocess.run(cmd, check=True)
        return subprocess.CompletedProcess(args=cmd, returncode=0)

    def run_command(self, tool: str, sample=None) -> list:
        """Generate commands for different tools"""
        BASE_PATH = "/Users/giuse/pythonProject/Mycodes/MyProjects"
        INPUT_PATH = f"{BASE_PATH}/simulated_ONT_data"
        RESULTS_PATH = f"{BASE_PATH}/results"
        CONDA_ENV = "long_reads" # name of your conda environment

        match tool:
            case "guppy_barcoder":
                cmd = [
                    "/Users/giuse/Tools/ont-guppy-cpu/bin/guppy_barcoder", # path of guppy tool
                    "-i", INPUT_PATH,
                    "-s", RESULTS_PATH,
                    "--trim_barcodes"
                ]

            case "porechop":
                cmd = ["conda", "run", "-n", CONDA_ENV, "porechop", "-i", INPUT_PATH]

            case "filtlong":
                input_file = f"{BASE_PATH}/data_separated/ONT_sample{sample}.fastq"
                output_file = f"{RESULTS_PATH}/barcode_{sample}.filtered.fastq"
                cmd = ["conda", "run",
                    "filtlong", "--min_length", "1000", "--keep_percent", "90",
                    "--trim", "-a", self.reference_genome, input_file, ">", output_file
                ]

            case "fastqc":
                input_file = f"{RESULTS_PATH}/barcode_{sample}.filtered.fastq"
                cmd = ["conda", "run", "-n", CONDA_ENV, "fastqc", input_file, "--outdir", f"{RESULTS_PATH}/fastqc"]

            case "Flye":
                cmd = ["conda", "run", "-n", CONDA_ENV,
                    "flye", "--nano-hq", f"{RESULTS_PATH}/barcode_{sample}.filtered.fastq",
                    "--genome-size", "4.6m", "--out-dir",
                    f"{RESULTS_PATH}/flye_assembly", "--threads", "7"
                ]

            case "Minimap2":
                input_file = f"{RESULTS_PATH}/barcode_{sample}.filtered.fastq"
                sam_output = f"{RESULTS_PATH}/sample_{sample}.sam"
                cmd = ["conda", "run", "-n", CONDA_ENV, "minimap2", "-ax", "map-ont", self.reference_genome, input_file, ">", sam_output]

            case "samtools_view":
                input_file = f"{RESULTS_PATH}/sample_{sample}.sam"
                bam_output = f"{RESULTS_PATH}/sample_{sample}.bam"
                cmd = ["conda", "run", "-n", CONDA_ENV, "samtools", "view", "-bS", input_file, "|", "samtools", "sort", "-o", bam_output]

            case "samtools_index":
                input_file = f"{RESULTS_PATH}/sample_{sample}.sam"
                bam_output = f"{RESULTS_PATH}/sample_{sample}.bam"
                cmd = ["conda", "run", "-n", CONDA_ENV, "samtools", "index", bam_output]

            case "picard":
                input_file = f"{RESULTS_PATH}/sample_{sample}.bam"
                output_file = f"{RESULTS_PATH}/sample_{sample}.markdup.bam"
                cmd = ["conda", "run", "-n", CONDA_ENV,
                    "picard", "MarkDuplicates", f"I={input_file}", f"O={output_file}",
                    f"METRICS_FILE={self.reference_genome}", "REMOVE_DUPLICATES=true", "CREATE_INDEX=true"
                ]

            case "CuteSV":
                input_file = f"{RESULTS_PATH}/sample_{sample}.markdup.bam"
                output_file = f"{RESULTS_PATH}/sample_{sample}.outcute"
                cmd = ["conda", "run", "-n", CONDA_ENV, "cuteSV", input_file, self.reference_genome, output_file]

            case "Sniffles":
                input_file = f"{RESULTS_PATH}/sample_{sample}.markdup.bam"
                output_file = f"{RESULTS_PATH}/sample_{sample}.vcf"
                cmd = ["conda", "run", "-n", CONDA_ENV, "sniffles", "--input", input_file, "--vcf", output_file]

            case _:
                raise ValueError(f"Unsupported tool: {tool}")

        return cmd


if __name__ == "__main__":
    BASE_PATH = "/Users/giuse/pythonProject/Mycodes/MyProjects" # set working directory
    FILE_PATH = f"{BASE_PATH}/simulated_ONT_data/combined.fastq" # file with all the reads
    SUMMARY_FILE_PATH = f"{BASE_PATH}/ONT_simulated_summary.txt" # summary output of ONT
    REFERENCE_GENOME_PATH = f"{BASE_PATH}/Escherichia_coli_reference" # reference genome

    for path in [FILE_PATH, SUMMARY_FILE_PATH, REFERENCE_GENOME_PATH]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"File not found: {path}")

    try:
        fastq_file = Fastq(FILE_PATH, REFERENCE_GENOME_PATH)
        run_quality_check(summary_file=SUMMARY_FILE_PATH)

        steps = [
            fastq_file.demultiplex(),
            fastq_file.filter_trim_reads(),
            fastq_file.reads_quality_check(),
            fastq_file.genome_assembly(),
            fastq_file.sort_index(),
            fastq_file.remove_duplicates(),
            fastq_file.variant_calling()
        ]

        for step in steps:
            if step.returncode != 0:
                print(f"Step failed: {step.args}")

    except subprocess.CalledProcessError as e:
        print(f"Tool failed: {e}")
    except ValueError as e:
        print(e)

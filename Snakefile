# Snakefile

# Configuration
configfile: "config.yaml"

rule all: # modify with the paths on your device
    input:
        "results/pycoqc_output.html",
        expand("results/guppy_output/barcode_{barcode}.fastq", barcode=range(1, 3)), 
        expand("results/filtered/barcode_{barcode}.filtered.fastq", barcode=range(1, 3)),
        expand("results/fastqc/barcode_{barcode}_report.html", barcode = range(1, 3)),
        expand("results/flye_assembly/assembly_{barcode}.fasta", barcode=range(1, 3)),
        expand("results/sample_{barcode}.bam", barcode=range(1, 3)),
        expand("results/sample_{barcode}.bai", barcode=range(1, 3)),
        expand("results/sample_{barcode}_marked.bam", barcode=range(1, 3)),
        expand("results/sample_{barcode}_cutesv.vcf", barcode=range(1, 3)),
        expand("results/sample_{barcode}_sniffles.vcf", barcode=range(1, 3))

rule pycoqc:
    input:
        config['summary_file']
    output:
        "results/pycoqc_output.html"
    shell:
        "conda run -n long_reads pycoQC --summary_file {input} -o {output}"

rule guppy_barcoder:
    input:
        config['input_fastq']
    output:
        expand("results/guppy_output/barcode_{barcode}.fastq", barcode=range(1, 3)) 
    shell:
        "/Users/giuse/Tools/ont-guppy-cpu/bin/guppy_barcoder -i {input} -s results/guppy_output" 

rule filtlong:
    input:
        "data_separated/barcode_{barcode}.fastq" # Example, adapt to your output
    output:
        "results/filtered/barcode_{barcode}.filtered.fastq"
    params:
        reference_genome = config['reference_genome']
    shell:
        "conda run filtlong --min_length 1000 --keep_percent 90 --trim -1 {params.reference_genome} {input} > {output}"

rule fastqc:
    input:
        "results/filtered/barcode_{barcode}.filtered.fastq"
    output:
        "results/fastqc/barcode_{barcode}_report.html" 
    shell:
        "conda run -n long_reads fastqc {input} --outdir results/fastqc/" 

rule flye:
    input:
        expand("results/filtered/barcode_{barcode}.filtered.fastq", barcode=range(1, 3))
    output:
        "results/flye_assembly/assembly_{barcode}.fasta"
    params:
        genome_size = config['genome_size']
    shell:
        "conda run -n long_reads flye --nano-hq {input} --genome-size {params.genome_size} --out-dir flye_assembly/ --threads 7"

rule minimap2:
    input:
        "results/filtered/barcode_{barcode}.filtered.fastq"
    output:
        "results/sample_{barcode}.sam"
    params:
        reference_genome = config['reference_genome']
    shell:
        "conda run -n long_reads minimap2 -ax map-ont {params.reference_genome} {input} > {output}"

rule samtools:
    input:
        "results/sample_{barcode}.sam"
    output:
        "results/sample_{barcode}.bam",
        "results/sample_{barcode}.bai"
    shell:
        """
        conda run -n long_reads samtools view -bS {input} | samtools sort -o {output[0]} 
        conda run -n long_reads samtools index {output[0]}
        """

rule picard_mark_duplicates:
    input:
        "results/sample_{barcode}.bam"
    output:
        "results/sample_{barcode}_marked.bam"
    shell:
        "conda run -n long_reads picard MarkDuplicates I={input} O={output} M=metrics.txt REMOVE_DUPLICATES=true" # Replace with your picard path

rule cutesv:
    input:
        "results/sample_{barcode}_marked.bam"
    output:
        "results/sample_{barcode}_cutesv.vcf"
    shell:
        "conda run -n long_reads cuteSV {input} {reference_genome} {output}"

rule sniffles:
    input:
        "results/sample_{barcode}_marked.bam"
    output:
        "results/sample_{barcode}_sniffles.vcf"
    shell:
        "conda run -n long_reads sniffles {input} {output}" 
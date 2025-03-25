import os.path
import vcfpy
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

def process_chunk(chunk):
    """Processes a chunk of VCF records into a DataFrame"""
    return pd.DataFrame(chunk)

class VCFProcessor:
    def __init__(self, input_vcfs : list, clinvar_vcf=None, gene_db=None):
        self.input_vcfs = input_vcfs if isinstance(input_vcfs, list) else [input_vcfs]
        self.merged_vcf = "merged_output.vcf"
        self.clinvar_vcf = clinvar_vcf
        self.gene_db = gene_db
        self.df = pd.DataFrame() # DataFrame to store data

    def merge_vcfs(self):
        """Merges multiple VCF into a single one"""
        if len(self.input_vcfs) == 1:
            print("Only one file found, skipping merging.")
            self.merged_vcf = self.input_vcfs[0]
            return

        print("Merging Files...")
        vcf_readers = [vcfpy.Reader(open(f, 'r')) for f in self.input_vcfs if os.path.exists(f)]
        if not vcf_readers:
            raise FileNotFoundError("No valid VCF files found.")
        vcf_writer = vcfpy.Writer(open(self.merged_vcf, 'w'), vcf_readers[0])

        for reader in vcf_readers:
            for record in reader:
                vcf_writer.write_record(record)

        print(f"Merged VCF saved to {self.merged_vcf}")

    def vcf_to_dataframe_parallel(self, num_workers=4):
        """Converts VCF into a Pandas DataFrame."""
        print("Converting VCF to DataFrame...")
        chunk_size = 50000
        vcf_reader = vcfpy.Reader(open(self.merged_vcf, 'r'))
        records = []

        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = []
            chunk = []

            for record in vcf_reader:
                chunk.append({
                    "CHROM": record.CHROM,
                    "POS": record.POS,
                    "ID": record.ID,
                    "REF":record.REF,
                    "ALT": ",".join(map(str, record.ALT)),
                    "QUAL": record.QUAL,
                    "FILTER": record.FILTER,
                    "DP": record.INFO.get("DP", 0),
                    "INFO": str(record.INFO)
                })

                if len(chunk) >= chunk_size:
                    futures.append(executor.submit(process_chunk, chunk))
                    chunk = [] # reset chunk

            if chunk:
                futures.append(executor.submit(process_chunk, chunk))

            df_chunks = [future.result() for future in as_completed(futures)]

        # Concatenate chunks
        self.df = pd.concat(df_chunks, ignore_index=True)
        print(f"Converted {len(self.df)} variants into DataFrame.")

    def filter_variants(self, min_qual=30, min_dp=10):
        """Filters variants based on QUAL and DP fields."""
        print("Filtering variants...")
        before_count = len(self.df)
        self.df = self.df[(self.df["QUAL"] >= min_qual) & (self.df["DP"] >= min_dp)]
        print(f"Remaining variants after filtering: {len(self.df)}")

    def count_variant_types(self):
        """Counts SNPs and Indels"""
        print("Counting SNPs and Indels...")
        self.df["Variant_Type"] = self.df.apply(lambda row: "SNP" if len(row["REF"]) == 1 and all(len(alt) == 1 for alt in row["ALT"].split(",")) else "INDEL", axis=1)
        counts = self.df["Variant_Type"].value_counts()
        print(counts)
        return counts

    def plot_qual_distribution(self):
        """Plots the distribution of variant quality scores."""
        print("Plotting variant quality distribution...")
        plt.hist(self.df["QUAL"], bins=50, color='blue', alpha=0.7)
        plt.xlabel("Variant Quality Score (QUAL)")
        plt.ylabel("Frequency")
        plt.title("Distribution of Variant Quality Scores")
        plt.show()

    def find_clinvar_variants(self):
        """Identifies variants present in ClinVar."""
        if not self.clinvar_vcf:
            print("No ClinVar file provided.")
            return

        print("Searching for ClinVar variants...")
        clinvar_df = pd.read_csv(self.clinvar_vcf, sep="\t")
        clinvar_positions = set(zip(clinvar_df["CHROM"], clinvar_df["POS"]))

        self.df["ClinVar_Match"] = self.df.apply(lambda row: (row["CHROM"], row["POS"]) in clinvar_positions, axis=1)
        matches = self.df[self.df["ClinVar_Match"]]
        print(f"Found {len(matches)} ClinVar matches.")
        return matches

    def annotate_variants(self):
        """Annotates variants file gene names."""
        if not self.gene_db:
            print("No gene database provided.")
            return

        print("Annotating variants...")
        gene_df = pd.read_csv(self.gene_db, sep="\t")
        gene_dict = dict(zip(zip(gene_df["CHROM"], gene_df["POS"]), gene_df["GENE"]))

        self.df["GENE"] = self.df.apply(lambda row:gene_dict.get((row["CHROM"], row["POS"]), "Unknown"), axis=1)
        print("Annotation completed.")

    def save_dataframe(self):
        """Save DataFrame to CSV"""
        self.df.to_csv("processed_variants.csv", index=False)
        print("Saved processed variants to CSV.")

    def run_pipeline(self):
        """Runs full VCF processing pipeline"""
        print("Starting VCF processing pipeline...")

        self.merge_vcfs()
        self.vcf_to_dataframe_parallel(num_workers=multiprocessing.cpu_count())
        self.filter_variants(min_qual=30, min_dp=10)
        self.count_variant_types()
        self.plot_qual_distribution()
        self.find_clinvar_variants()
        self.annotate_variants()
        self.save_dataframe()

        print("VCF processing completed.")

# Run Pipeline
if __name__ == "__main__":
    input_vcf = ["/Users/giuse/Downloads/case584_variants.vcf"]
    clinvar_vcf = "clinvar.vcf"
    gene_db = "gene_database.txt"

    processor = VCFProcessor(input_vcf, clinvar_vcf, gene_db)
    processor.run_pipeline()

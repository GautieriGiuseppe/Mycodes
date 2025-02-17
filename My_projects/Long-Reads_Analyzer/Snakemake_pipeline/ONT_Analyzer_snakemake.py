from snakemake import snakemake

# Define path of Snakefile
snakefile_path = ""

# run snakemake pipeline
success = snakemake(
    snakefile_path,
    cores=4,
    keepgoing=False, # Stop if something fails
    forceall=False, # Set True if want to rerun all rules
    printshellcmds=True # Print shell commands as they are executed
)

# Check if Snakemake completed successfully
if success:
    print("Snakemake pipeline completed successfully!")
else:
    print("Snakemake pipeline failed.")
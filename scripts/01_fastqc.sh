#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=1000M
#SBATCH --time=3:00:00
#SBATCH --job-name=fastqc
#SBATCH --mail-user=anmol.ratan@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/data/users/aratan/rnaseq/fastqc/output_fastqc_%j.o
#SBATCH --error=/data/users/aratan/rnaseq/fastqc/error_fastqc_%j.e
#SBATCH --partition=pibu_el8

# Define the new working directory containing FASTQ files
READS_DIR="/data/users/aratan/rnaseq/breastcancer_de/reads"

# Output directory for FastQC results
OUTPUT_DIR="/data/users/aratan/rnaseq/fastqc"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Run FastQC using Apptainer
apptainer exec \
--bind $READS_DIR:$READS_DIR --bind $OUTPUT_DIR:$OUTPUT_DIR \
/containers/apptainer/fastqc-0.12.1.sif \
bash -c "fastqc -o $OUTPUT_DIR $READS_DIR/*.fastq.gz"


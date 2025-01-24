#!/bin/bash
#SBATCH --job-name=sam_to_bam_array         # Job name
#SBATCH --partition=pibu_el8               # Partition name
#SBATCH --cpus-per-task=1                  # Number of CPU cores
#SBATCH --mem=4000                         # Memory allocation
#SBATCH --time=1:00:00                     # Maximum runtime
#SBATCH --mail-user=anmol.ratan@students.unibe.ch
#SBATCH --mail-type=END
#SBATCH --array=0-11                       # Array index range (adjust based on number of samples)
#SBATCH --output=/data/users/aratan/rnaseq/map/bam/output_sam_to_bam_%A_%a.o
#SBATCH --error=/data/users/aratan/rnaseq/map/bam/error_sam_to_bam_%A_%a.e

# Directories
SAM_DIR="/data/users/aratan/rnaseq/map/mapped_reads"
BAM_DIR="/data/users/aratan/rnaseq/map/bam"
CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

mkdir -p $BAM_DIR

SAMPLES=(HER21 HER22 HER23 TNBC1 TNBC2 TNBC3 NonTNBC1 NonTNBC2 NonTNBC3 Normal1 Normal2 Normal3)

# Get the current sample
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
SAM_FILE=${SAM_DIR}/${SAMPLE}.sam
BAM_FILE=${BAM_DIR}/${SAMPLE}.bam

# Convert SAM to BAM using the container
apptainer exec --bind /data/ $CONTAINER \
samtools view -hbS $SAM_FILE > $BAM_FILE

echo "Converted $SAM_FILE to $BAM_FILE"


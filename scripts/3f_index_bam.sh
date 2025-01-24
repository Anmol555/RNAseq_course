#!/bin/bash
#SBATCH --job-name=index_all_bams         # Job name
#SBATCH --output=/data/users/aratan/rnaseq/map/indexbam/output_index_bam_%A_%a.o  # Standard output
#SBATCH --error=/data/users/aratan/rnaseq/map/indexbam/error_index_bam_%A_%a.e   # Standard error
#SBATCH --mail-user=anmol.ratan@students.unibe.ch  # Email for notifications
#SBATCH --mail-type=END                          # Send email at job completion
#SBATCH --time=01:00:00                          # Maximum runtime (1 hour)
#SBATCH --ntasks=1                               # Number of tasks (1 CPU)
#SBATCH --cpus-per-task=1                        # Number of CPU cores per task
#SBATCH --mem=4000                               # Memory in MB (4 GB)
#SBATCH --partition=pibu_el8                        # Partition (adjust as needed)
#SBATCH --array=0-11                             # Array job for 12 files (update this as per your BAM file count)

# Define variables
CONTAINER_PATH="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
DATA_DIR="/data/users/aratan/rnaseq/map/sorted_bam/"

# Get the BAM file name based on SLURM_ARRAY_TASK_ID
BAM_FILES=($(ls $DATA_DIR/*.bam))
BAM_FILE=${BAM_FILES[$SLURM_ARRAY_TASK_ID]}

# Run samtools index using Apptainer
apptainer exec --bind $DATA_DIR $CONTAINER_PATH samtools index $BAM_FILE


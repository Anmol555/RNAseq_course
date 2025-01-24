#!/bin/bash
#SBATCH --job-name=sort_bam_array         # Job name
#SBATCH --partition=pibu_el8             # Partition name
#SBATCH --cpus-per-task=4                # 4 CPUs per task
#SBATCH --mem=64G                      # Memory allocation (25,000 MB)
#SBATCH --time=10:00:00                   # Maximum runtime (2 hours)
#SBATCH --array=0-11                     # Array index range (adjust based on the number of BAM files)
#SBATCH --mail-user=anmol.ratan@students.unibe.ch
#SBATCH --mail-type=END
#SBATCH --output=/data/users/aratan/rnaseq/map/sorted_bam/output_sort_bam_%A_%a.o
#SBATCH --error=/data/users/aratan/rnaseq/map/sorted_bam/error_sort_bam_%A_%a.e

# Define directories
INPUT_DIR="/data/users/aratan/rnaseq/map/bam"
OUTPUT_DIR="/data/users/aratan/rnaseq/map/sorted_bam"
TEMP_DIR="/data/users/aratan/rnaseq/map/bam/temp"

# Create temporary and output directories if they don't exist
mkdir -p $TEMP_DIR
mkdir -p $OUTPUT_DIR

# Get the list of BAM files dynamically
FILES=($(ls $INPUT_DIR/*.bam))

# Determine the current BAM file based on SLURM_ARRAY_TASK_ID
CURRENT_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
BASENAME=$(basename $CURRENT_FILE .bam)
OUTPUT_FILE="$OUTPUT_DIR/${BASENAME}_sorted.bam"

# Run Apptainer with Samtools for sorting
echo "Sorting file: $CURRENT_FILE"
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif samtools sort \
    -m 6G -@ 4 -o $OUTPUT_FILE -T $TEMP_DIR/$BASENAME $CURRENT_FILE

# Completion message
if [[ $? -eq 0 ]]; then
    echo "Sorting of $CURRENT_FILE completed successfully. Output saved as $OUTPUT_FILE"
else
    echo "Error occurred while sorting $CURRENT_FILE"
    exit 1
fi


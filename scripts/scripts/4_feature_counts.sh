#!/bin/bash

#SBATCH --job-name=featureCounts_job         # Job name
#SBATCH --partition=pibu_el8                # Replace with your cluster's partition
#SBATCH --ntasks=1                          # Number of tasks
#SBATCH --cpus-per-task=4                   # CPU cores (as required)
#SBATCH --mem=4000MB                        # Memory (4GB as required)
#SBATCH --time=13:00:00                     # Time limit (13 hours)
#SBATCH --mail-user=anmol.ratan@students.unibe.ch  # Email address for notifications
#SBATCH --mail-type=END                     # Send email when the job ends
#SBATCH --output=/data/users/aratan/rnaseq/count/output_featureCounts_%A_%a.o
#SBATCH --error=/data/users/aratan/rnaseq/count/error_featureCounts_%A_%a.e


# Define paths
CONTAINER="/containers/apptainer/subread_2.0.1--hed695b0_0.sif"  # Path to Apptainer container
ANNOTATION="/data/users/aratan/rnaseq/map/Homo_sapiens.GRCh38.113.gtf"  # Annotation file
BAM_DIR="/data/users/aratan/rnaseq/map/sorted_bam"   # Directory containing sorted BAM files
OUTPUT="/data/users/aratan/rnaseq/count/gene_counts.txt"  # Output file for counts


# Decompress GTF file if it is in .gz format
if [[ -f "${ANNOTATION}.gz" ]]; then
    echo "Decompressing annotation file..."
    gunzip -f "${ANNOTATION}.gz"
fi


# Run featureCounts with bind mounts
apptainer exec --bind /data:/data $CONTAINER featureCounts -T 4 -Q 10 -a $ANNOTATION -o $OUTPUT $BAM_DIR/*.bam

# Print completion message
echo "FeatureCounts job completed successfully!"

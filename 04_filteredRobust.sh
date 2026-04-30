#!/bin/bash

#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=filter_ITS
#SBATCH --mail-type=begin,end
#SBATCH --time=04:00:00
#SBATCH --mail-user=adamantia.kapopoulou@unibe.ch
#SBATCH --output=seqkit_%j.out
#SBATCH --error=seqkit_%j.err
#SBATCH --partition=pibu_el8

module load SeqKit/2.6.1

# Threads (safe fallback)
THREADS=${SLURM_CPUS_PER_TASK:-4}

# Input/output directories
INPUT_DIR=/data/projects/p605_PerchPilot/PacBio/Colin/ITS/raw/
OUTPUT_DIR=/data/projects/p605_PerchPilot/PacBio/Colin/ITS/raw/filteredRobust/

seqkit stats ${INPUT_DIR}/*.fastq.gz > ${OUTPUT_DIR}/stats_before.txt

# Loop through all fastq.gz files
for file in ${INPUT_DIR}/*.fastq.gz
do
	base=$(basename ${file} .fastq.gz)
	echo "[$(date)] Filtering: $base"
	seqkit seq ${file} -M 2000 -m 200 -j ${THREADS} |
	seqkit grep -s -v -p "N" -j ${THREADS} -o ${OUTPUT_DIR}/${base}.filtered.fastq.gz
done
echo "All samples processed."

seqkit stats ${OUTPUT_DIR}/*.filtered.fastq.gz > ${OUTPUT_DIR}/stats_after.txt
echo "Step 2 done."

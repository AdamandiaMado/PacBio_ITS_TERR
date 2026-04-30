#!/bin/bash

#SBATCH --mem-per-cpu=16G
#SBATCH --ntasks=8
#SBATCH --job-name=seqkit_ITS
#SBATCH --mail-type=begin,end
#SBATCH --time=04:00:00
#SBATCH --mail-user=adamantia.kapopoulou@unibe.ch
#SBATCH --output=seqkit_%j.out
#SBATCH --error=seqkit_%j.err
#SBATCH --partition=pibu_el8

module load SeqKit/2.6.1

# Input/output directories
INPUT_DIR=/data/projects/p605_PerchPilot/PacBio/Colin/ITS/raw
OUTPUT_DIR=/data/projects/p605_PerchPilot/PacBio/Colin/ITS/raw/filtered/
# Loop through all fastq.gz files
for file in ${INPUT_DIR}/*.fastq.gz
do
	base=$(basename ${file} .fastq.gz)
	echo "Processing ${base}..."
	seqkit stats ${file} -o ${OUTPUT_DIR}/${base}
done
echo "All samples processed."

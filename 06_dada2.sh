#!/bin/bash 

#SBATCH --mem-per-cpu=54G 
#SBATCH --cpus-per-task=6
#SBATCH --job-name=dada2_ITS 
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=adamantia.kapopoulou@unibe.ch
#SBATCH --time=69:00:00 
#SBATCH --output=dada2_%j.out 
#SBATCH --error=dada2_%j.err 
#SBATCH --partition=pibu_el8 

module load R

INPUT_DIR=/data/projects/p605_PerchPilot/PacBio/Colin/ITS/raw/filteredRobust/downsampled
OUTPUT_DIR=/data/projects/p605_PerchPilot/PacBio/Colin/ITS/raw/filteredRobust/downsampled/dada2/ 
SCRIPT=06b_resumeDada2.R 

Rscript ${SCRIPT} ${INPUT_DIR} ${OUTPUT_DIR}

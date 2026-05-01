#!/bin/bash
#SBATCH --mem-per-cpu=32G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=taxonomy
#SBATCH --mail-type=begin,end
#SBATCH --time=28:00:00
#SBATCH --output=taxonomy_%j.out
#SBATCH --error=taxonomy_%j.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=adamantia.kapopoulou@unibe.ch

module load R-bundle-Bioconductor/3.18-foss-2021a-R-4.3.2
unset R_LIBS_USER

Rscript taxonomy_phyloseq.R

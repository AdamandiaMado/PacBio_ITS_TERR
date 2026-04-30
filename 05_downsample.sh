#!/bin/bash
#SBATCH --job-name=downsample_ITS
#SBATCH --output=logs/downsample_%A_%a.out
#SBATCH --error=logs/downsample_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-type=begin,end
#SBATCH --partition=pibu_el8
#SBATCH --array=1-211
#SBATCH --mail-user=adamantia.kapopoulou@unibe.ch

module load SeqKit/2.6.1

INPUT_DIR=/data/projects/p605_PerchPilot/PacBio/Colin/ITS/raw/filteredRobust
OUTPUT_DIR=/data/projects/p605_PerchPilot/PacBio/Colin/ITS/raw/filteredRobust/downsampled
FILES=(${INPUT_DIR}/*filtered.fastq.gz)
f=${FILES[$SLURM_ARRAY_TASK_ID-1]}
base=$(basename "$f" .filtered.fastq.gz)
echo "Downsampling $base (task ${SLURM_ARRAY_TASK_ID})..."
N_READS=$(seqkit stats "$f" | awk 'NR==2{print $4}' | tr -d ',')

if [[ "${N_READS}" -lt 20000 ]]; then
    echo "  WARNING: only ${N_READS} reads — copying as-is (below 20,000 threshold)"
    cp "$f" "${OUTPUT_DIR}/${base}.downsampled.fastq.gz"
else
    echo "  ${N_READS} reads → downsampling to 20,000"
    seqkit sample -n 20000 -s 42 "$f" -o "${OUTPUT_DIR}/${base}.downsampled.fastq.gz"
fi
echo "Done $base"

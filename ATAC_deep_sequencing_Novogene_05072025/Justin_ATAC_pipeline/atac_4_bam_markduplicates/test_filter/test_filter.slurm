#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=384G
#SBATCH --time=24:00:00
#SBATCH --verbose
#SBATCH -J picard.markdups
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# example of sbatch command to run this script:
# sbatch -p hpc_l40s --array=6,8,9 --export=inputfiles=../inputfiles.txt,outdir=$(pwd),samplelist=../samplenames.txt test_filter.slurm

INPUT=$(head -n $SLURM_ARRAY_TASK_ID "$inputfiles" | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID "$samplelist" | tail -n 1)
OUTPUT_DIR="${outdir}/${NAME}"

mkdir -p "$OUTPUT_DIR"

source activate fastq2bam

set -euo pipefail
# Step 1: Filter to keep properly paired reads. -F 4 to discard only unmapped reads but keeps duplicates and reads with secondary alignment since we are doing centromeric analysis. if I add the usual -q 30 (mapping quality threshold) it will throw out most of the reads in this case, since multi-mappers usually have q=0, so exclude the quality thresholding here. Sort for fixmate
samtools view -@ 6 -F 4 -f 2 -u "$INPUT" | \
  samtools sort -n --threads 6 - > "${OUTPUT_DIR}/${NAME}_f2.bam"
  
samtools flagstat -@ 6 "${OUTPUT_DIR}/${NAME}_f2.bam" > "${OUTPUT_DIR}/${NAME}_f2.flagstat"

samtools view -@ 6 -F 4 -f 1 -u "$INPUT" | \
  samtools sort -n --threads 6 - > "${OUTPUT_DIR}/${NAME}_f1.bam"
  
samtools flagstat -@ 6 "${OUTPUT_DIR}/${NAME}_f1.bam" > "${OUTPUT_DIR}/${NAME}_f1.flagstat"

samtools view -@ 6 -F 4 -u "$INPUT" | \
  samtools sort -n --threads 6 - > "${OUTPUT_DIR}/${NAME}_no_f.bam"
  
samtools flagstat -@ 6 "${OUTPUT_DIR}/${NAME}_no_f.bam" > "${OUTPUT_DIR}/${NAME}_no_f.flagstat"


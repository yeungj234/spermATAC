#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=384G
#SBATCH --time=24:00:00
#SBATCH --verbose
#SBATCH -J picard.markdups
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# example of sbatch command to run this script:
# sbatch -p hpc_l40_b --array=1-6 --export=inputfiles=inputfiles.txt,outdir=$(pwd),samplelist=samplenames.txt atac_4_bam_markduplicates.slurm

INPUT=$(head -n $SLURM_ARRAY_TASK_ID "$inputfiles" | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID "$samplelist" | tail -n 1)
OUTPUT_DIR="${outdir}/${NAME}"

mkdir -p "$OUTPUT_DIR"

source activate fastq2bam

set -euo pipefail
# Step 1: Filter to keep properly paired reads. -F 4 to discard only unmapped reads but keeps duplicates and reads with secondary alignment since we are doing centromeric analysis. if I add the usual -q 30 (mapping quality threshold) it will throw out most of the reads in this case, since multi-mappers usually have q=0, so exclude the quality thresholding here. Sort for fixmate
#samtools view -@ 16 -F 4 -f 2 -u "$INPUT" | \
#  samtools sort -n --threads 16 - > "${OUTPUT_DIR}/${NAME}_filter.bam"
  
#samtools flagstat -@ 16 "${OUTPUT_DIR}/${NAME}_filter.bam" > "${OUTPUT_DIR}/${NAME}_filter.flagstat"

# Step 2: Fixmate: Adds mate information (like insert size, mate position) and flags proper pairs. Removed the -r argument since we want to keep secondary alignments since we are doing centromeric analysis. 
#samtools fixmate -@ 16 "${OUTPUT_DIR}/${NAME}_filter.bam" "${OUTPUT_DIR}/${NAME}_fixmate.bam"
#samtools flagstat -@ 16 "${OUTPUT_DIR}/${NAME}_fixmate.bam" > "${OUTPUT_DIR}/${NAME}_fixmate.flagstat"
#samtools view -@ 16 -u "${OUTPUT_DIR}/${NAME}_fixmate.bam" | samtools sort --threads 16 - -o "${OUTPUT_DIR}/${NAME}_fixmate_sort.bam"

# Cleanup
# rm "${OUTPUT_DIR}/${NAME}_filter.bam"
 rm "${OUTPUT_DIR}/${NAME}_fixmate.bam"
 
# Step 3: Mark duplicates
picard MarkDuplicates \
  I="${OUTPUT_DIR}/${NAME}_fixmate_sort.bam" \
  O="${OUTPUT_DIR}/${NAME}_markdups.bam" \
  M="${OUTPUT_DIR}/${NAME}_metrics.txt"

# Step 4: Sort and index final bam
samtools sort -@ 16 "${OUTPUT_DIR}/${NAME}_markdups.bam" -o "${OUTPUT_DIR}/${NAME}_markdups_final.bam"
samtools index -@ 16 "${OUTPUT_DIR}/${NAME}_markdups_final.bam"
samtools flagstat -@ 16 "${OUTPUT_DIR}/${NAME}_markdups_final.bam" > "${OUTPUT_DIR}/${NAME}_markdups_final.flagstat"

# Step 5: split the final bam file by overlapping vs not overlapping with centromere (see if samtools view -L option doesn't run out of memory) 
# Get overlapping reads
samtools view -@ 16 -b -L /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/peaks/GSE153058_xla_v10.2_cen.bed "${OUTPUT_DIR}/${NAME}_markdups_final.bam" > ${OUTPUT_DIR}/${NAME}_centromere.bam
# Get non-overlapping reads & discard overlapping reads. 
samtools view -@ 16 -b -L /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/peaks/GSE153058_xla_v10.2_cen.bed -U "${OUTPUT_DIR}/${NAME}_noncentromere.bam" "${OUTPUT_DIR}/${NAME}_markdups_final.bam" > /dev/null

# Step 6: Compute Centromeric Read Coverage
bedtools coverage -sorted \
-a /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/peaks/GSE153058_xla_v10.2_cen.bed \
-b "${OUTPUT_DIR}/${NAME}_markdups_final.bam" > ${OUTPUT_DIR}/${NAME}_centromere_coverage.txt -hist

# Step 6: split the final bam file by overlapping vs not overlapping with centromere
# bedtools intersect -abam "${OUTPUT_DIR}/${NAME}_markdups_final.bam" -b /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/peaks/GSE153058_xla_v10.2_cen.bed -ubam -wa > ${OUTPUT_DIR}/${NAME}_centromere.bam
# bedtools subtract -abam "${OUTPUT_DIR}/${NAME}_markdups_final.bam" -b /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/peaks/GSE153058_xla_v10.2_cen.bed -ubam -wa > ${OUTPUT_DIR}/${NAME}_noncentromere.bam

# Step 7: 
# Cleanup
 rm "${OUTPUT_DIR}/${NAME}_fixmate_sort.bam"
 rm "${OUTPUT_DIR}/${NAME}_markdups.bam"
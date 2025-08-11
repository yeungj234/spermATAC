#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=384G
#SBATCH --time=24:00:00
#SBATCH --verbose
#SBATCH -J centromere.ol
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# Example of sbatch command used
# sbatch -p hpc_l40s_b --array=1-6 --export=outdir=$(pwd),samplelist=samplenames.txt bedtools_centromeric_overlap.slurm

NAME=$(head -n $SLURM_ARRAY_TASK_ID "$samplelist" | tail -n 1)
OUTPUT_DIR="${outdir}/${NAME}"

source activate fastq2bam

set -euo pipefail

# Step 4: Sort and index final bam
samtools sort -@ 6 "${OUTPUT_DIR}/${NAME}_markdups.bam" -o "${OUTPUT_DIR}/${NAME}_markdups_final.bam"
samtools index -@ 6 "${OUTPUT_DIR}/${NAME}_markdups_final.bam"
samtools flagstat -@ 6 "${OUTPUT_DIR}/${NAME}_markdups_final.bam" > "${OUTPUT_DIR}/${NAME}_markdups_final.flagstat"

# Step 5: split the final bam file by overlapping vs not overlapping with centromere
# bedtools intersect -abam "${OUTPUT_DIR}/${NAME}_markdups_final.bam" -b /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/peaks/GSE153058_xla_v10.2_cen.bed -ubam -wa > ${OUTPUT_DIR}/${NAME}_centromere.bam
# bedtools subtract -abam "${OUTPUT_DIR}/${NAME}_markdups_final.bam" -b /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/peaks/GSE153058_xla_v10.2_cen.bed -ubam -wa > ${OUTPUT_DIR}/${NAME}_noncentromere.bam

# Step 5: split the final bam file by overlapping vs not overlapping with centromere (see if samtools view -L option doesn't run out of memory) 
# Get overlapping reads
samtools view -b -L /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/peaks/GSE153058_xla_v10.2_cen.bed "${OUTPUT_DIR}/${NAME}_markdups_final.bam" > ${OUTPUT_DIR}/${NAME}_centromere.bam
# Get non-overlapping reads & discard overlapping reads. 
samtools view -b -L /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/peaks/GSE153058_xla_v10.2_cen.bed -U "${OUTPUT_DIR}/${NAME}_noncentromere.bam" "${OUTPUT_DIR}/${NAME}_markdups_final.bam" > /dev/null

# Step 6: Compute Centromeric Read Coverage
bedtools coverage -sorted \
-a /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/peaks/GSE153058_xla_v10.2_cen.bed \
-b "${OUTPUT_DIR}/${NAME}_markdups_final.bam" > ${OUTPUT_DIR}/${NAME}_centromere_coverage.txt -hist

# Step 7: 
# Cleanup
 rm "${OUTPUT_DIR}/${NAME}_fixmate_sort.bam"
 rm "${OUTPUT_DIR}/${NAME}_markdups_final.bam"

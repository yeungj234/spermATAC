#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --verbose
#SBATCH -J SortIndex
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# sbatch -p hpc_l40_b --array=1-2 --export=inputfiles=sort_input.txt,outdir=$(pwd),samplelist=sort_samplenames.txt sort_index.slurm
# Read input and sample name
INPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$inputfiles")
NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$samplelist")
OUTPUT_DIR="${outdir}/${NAME}"

source activate fastq2bam

# samtools sort -@ 4 -o ./atacHAK034_WT/atacHAK034_WT_sorted.bam ./atacHAK034_WT/atacHAK034_WT_1447257.bam
# samtools sort -@ 4 -o ./atacHAK035_43AMut/atacHAK035_sorted.bam ./atacHAK035_43AMut/atacHAK035_43AMut_1447258.bam
# samtools sort -@ 4 -o ./atacHAK036_HIDMut/atacHAK036_sorted.bam ./atacHAK036_HIDMut/atacHAK036_HIDMut_1447256.bam
# samtools sort -@ 4 -o ./atacHAK040_WT/atacHAK040_WT_sorted.bam ./atacHAK040_WT/atacHAK040_WT_1447321.bam
# samtools sort -@ 4 -o ./atacHAK041_43AMut/atacHAK041_43AMut_sorted.bam ./atacHAK041_43AMut/atacHAK041_43AMut_1447322.bam
# samtools sort -@ 4 -o ./atacHAK042_HIDMut/atacHAK042_HIDMut_sorted.bam ./atacHAK042_HIDMut/atacHAK042_HIDMut_1447802.bam
# samtools sort -@ 4 -o ./atacHAK043_WT/atacHAK043_WT_sorted.bam ./atacHAK043_WT/atacHAK043_WT_1448447.bam

samtools sort -@ 8 -o ${OUTPUT_DIR}/${NAME}_sorted.bam $INPUT
# samtools index ./atacHAK034_WT/atacHAK034_WT_sorted.bam
# samtools index ./atacHAK035_43AMut/atacHAK035_sorted.bam
# samtools index ./atacHAK036_HIDMut/atacHAK036_sorted.bam
# samtools index ./atacHAK040_WT/atacHAK040_WT_sorted.bam
# samtools index ./atacHAK041_43AMut/atacHAK041_43AMut_sorted.bam
# samtools index ./atacHAK042_HIDMut/atacHAK042_HIDMut_sorted.bam
# samtools index ./atacHAK043_WT/atacHAK043_WT_sorted.bam

samtools index ${OUTPUT_DIR}/${NAME}_sorted.bam
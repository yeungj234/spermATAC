#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --verbose
#SBATCH -J manageMito
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# Example sbatch command:
# sbatch -p hpc_l40s_b --array=1-6 --export=inputfiles=inputfiles.txt,outdir=$(pwd),samplelist=samplenames.txt atac_3_bam_managemito.slurm

# Read input and sample name
INPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$inputfiles")
NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$samplelist")
OUTPUT_DIR="${outdir}/${NAME}"
mkdir -p "$OUTPUT_DIR"

source /ru-auth/local/home/risc_soft/miniconda3/etc/profile.d/conda.sh
conda activate /ru-auth/local/home/risc_soft/miniconda3/envs/fastq2bam

set -euo pipefail

# Input contig lists
chromNotMito="/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Justin_ATAC_pipeline/atac_3_bam_managemito/chromNotMito.txt"
chrom122XY="/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Justin_ATAC_pipeline/atac_3_bam_managemito/chrom1-10.txt"

# Temp dir for samtools sort
TMPDIR="${OUTPUT_DIR}/samtools_tmp"
mkdir -p "$TMPDIR"

# Extract non-mito BAM
samtools view -@ 4 -b "$INPUT" $(cat "$chromNotMito") > "${OUTPUT_DIR}/${NAME}_nonMito.bam"

# Extract chr1-10 BAM
samtools view -@ 4 -b "$INPUT" $(cat "$chrom122XY") > "${OUTPUT_DIR}/${NAME}_chr1-10.bam"

# Flagstat for nonMito BAM
samtools sort -n --threads 4 -T "${TMPDIR}/nonMito_sort" "${OUTPUT_DIR}/${NAME}_nonMito.bam" -O SAM | \
  samtools flagstat -@ 4 - > "${OUTPUT_DIR}/${NAME}_nonMito.flagstat"

# Flagstat for chr1-10 BAM
samtools sort -n --threads 4 -T "${TMPDIR}/chr_sort" "${OUTPUT_DIR}/${NAME}_chr1-10.bam" -O SAM | \
  samtools flagstat -@ 4 - > "${OUTPUT_DIR}/${NAME}_chr1-10.flagstat"

# Clean up temp dir
rm -rf "$TMPDIR"
rm "${OUTPUT_DIR}/${NAME}_nonMito.bam"

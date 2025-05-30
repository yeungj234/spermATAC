#!/bin/bash
#SBATCH --job-name=kmer_array
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --verbose
#SBATCH --cpus-per-task=16

# example of sbatch command
# sbatch --array=1-6 --export=inputfiles=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Justin_ATAC_pipeline/atac_2_fastq_align/inputfiles_R1.txt,samplelist=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Justin_ATAC_pipeline/atac_2_fastq_align/samplenames.txt,outdir=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Centromeric_Repeat_Seq/count_kmers_in_fastq_verbose/output/10mers array_10mer_count.sh

# Initialize conda properly
source ~/.bashrc
conda activate /lustre/fs4/risc_lab/store/jyeung/miniconda3/envs/pyahocorasick

# Get the input file and sample name for this array task
INPUT=$(head -n $SLURM_ARRAY_TASK_ID "$inputfiles" | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID "$samplelist" | tail -n 1)

# Set output file paths
OUT_COUNT="${outdir}/${NAME}_counts.tsv"
OUT_DETAIL="${outdir}/${NAME}_matches.tsv"

# Run the k-mer matching Python script
python /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Centromeric_Repeat_Seq/count_kmers_in_fastq_verbose/count_kmers_in_fastq_verbose.py \
  -k /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Centromeric_Repeat_Seq/FCR_monomers/FCR_10mers_seqonly.fasta \
  -f "$INPUT" \
  -o "$OUT_COUNT" \
  -d "$OUT_DETAIL"

#!/bin/bash
#SBATCH --ntasks=2                  
#SBATCH --cpus-per-task=16        
#SBATCH --time=24:00:00  
#SBATCH --mem=64G       
#SBATCH --verbose
#SBATCH -J clumpify
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# example of sbatch command to run this script:
# sbatch --array=1-6 --export=inputfiles_R1=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Justin_ATAC_pipeline/atac_2_fastq_align/inputfiles_R1.txt,inputfiles_R2=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Justin_ATAC_pipeline/atac_2_fastq_align/inputfiles_R2.txt,outdir=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Centromeric_Repeat_Seq/clumpify,samplelist=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Justin_ATAC_pipeline/atac_2_fastq_align/samplenames.txt clumpify_dedup.slurm

READ1=$(head -n $SLURM_ARRAY_TASK_ID "$inputfiles_R1" | tail -n 1)
READ2=$(head -n $SLURM_ARRAY_TASK_ID "$inputfiles_R2" | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID "$samplelist" | tail -n 1)
OUTDIR=$outdir\/$NAME

source activate Hide_ATAC
cd $outdir
# Set a temporary working file location
TMP_OUT1="${OUTDIR}/${NAME}_R1.dedup.fastq.gz"
TMP_OUT2="${OUTDIR}/${NAME}_R2.dedup.fastq.gz"

# Run Clumpify to remove PCR duplicates before trying to find kmers from FCR monomers are found in the fastqs. 
clumpify.sh \
  in1="$READ1" in2="$READ2" \
  out1="$TMP_OUT1" out2="$TMP_OUT2" \
  dedupe > /dev/null
  
# Run seqkit stats
seqkit stats "$READ1" > "${OUTDIR}/${NAME}_seqkit_stats.txt"
seqkit stats "$TMP_OUT1" >> "${OUTDIR}/${NAME}_seqkit_stats.txt"
seqkit stats "$READ2" >> "${OUTDIR}/${NAME}_seqkit_stats.txt"
seqkit stats "$TMP_OUT2" >> "${OUTDIR}/${NAME}_seqkit_stats.txt"

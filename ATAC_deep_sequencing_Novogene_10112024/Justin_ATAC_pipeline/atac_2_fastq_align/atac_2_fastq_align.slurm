#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --verbose
#SBATCH -J bowtie2
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# example of sbatch command to run this script:
# input fastqs are fastqs after trimming off adapters and PCR duplicate removal using clumpify. 
# sbatch -p hpc_a10_a --array=1-6 --export=inputfiles_R1=inputfiles_R1.txt,inputfiles_R2=inputfiles_R2.txt,outdir=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Justin_ATAC_pipeline/atac_2_fastq_align_05272025,genomeIndexName=/rugpfs/fs0/risc_lab/store/risc_data/downloaded/Xenla10.1/genome/Sequence/Bowtie2Index/genome,samplelist=samplenames.txt atac_2_fastq_align.slurm
INPUT_R1=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles_R1 | tail -n 1)
INPUT_R2=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles_R2 | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID $samplelist | tail -n 1)
OUTPUT_DIR=${outdir}/${NAME}
mkdir -p "$OUTPUT_DIR"

source activate fastq2bam

# Output BAM filename
BAM_OUT="${OUTPUT_DIR}/${NAME}_${SLURM_JOB_ID}.bam"

# Run Bowtie2 -> samtools view -> samtools sort
bowtie2 -k 100 -p 16 -X2000 --local --mm -x "$genomeIndexName" \
  -1 "$INPUT_R1" -2 "$INPUT_R2" --seed 23 2> "${BAM_OUT%.bam}.log" \
|samtools view -Su -@ 16 - \
|samtools sort -n -@ 16 -o "$BAM_OUT" -

# Index the BAM file
samtools index -@ 16 "$BAM_OUT"

# Run flagstat to generate QC stats
samtools flagstat -@ 16 "$BAM_OUT" > "${BAM_OUT%.bam}_flagstat.txt"

#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=384G
#SBATCH --time=24:00:00
#SBATCH --verbose
#SBATCH -J bowtie2
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# example of sbatch command to run this script:
# input fastqs are fastqs after trimming off adapters and PCR duplicate removal using clumpify. 
# sbatch --dependency=afterok:1446041 -p hpc_a10_a --array=1-9 --export=inputfiles_R1=inputfiles_R1.txt,inputfiles_R2=inputfiles_R2.txt,outdir=$(pwd),genomeIndexName=/rugpfs/fs0/risc_lab/store/risc_data/downloaded/Xenla10.1/genome/Sequence/Bowtie2Index/genome,samplelist=samplenames.txt atac_2_fastq_align.slurm

INPUT_R1=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles_R1 | tail -n 1)
INPUT_R2=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles_R2 | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID $samplelist | tail -n 1)
OUTPUT_DIR=${outdir}/${NAME}
mkdir -p "$OUTPUT_DIR"

source activate fastq2bam

# Output BAM filename
BAM_OUT="${OUTPUT_DIR}/${NAME}_${SLURM_JOB_ID}.bam"

# Run Bowtie2 -> samtools view -> samtools sort
bowtie2 -k 20 -p 16 -X2000 --local --mm -x "$genomeIndexName" \
  -1 "$INPUT_R1" -2 "$INPUT_R2" --seed 23 2> "${BAM_OUT%.bam}.log" \
|samtools view -Su -@ 8 - \
|samtools sort -n -@ 8 -m 4G -o "$BAM_OUT" -

# Index the BAM file
samtools index "$BAM_OUT"

# Run flagstat to generate QC stats
samtools flagstat -@ 8 "$BAM_OUT" > "${BAM_OUT%.bam}_flagstat.txt"

#! /bin/bash
#SBATCH --verbose
#SBATCH -J bowtie2_EColi
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out
source activate fastq2bam
INPUT_R1=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles_R1 | tail -n 1)
INPUT_R2=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles_R2 | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID $samplelist | tail -n 1)
OUTPUT=$outdir\/bowtie_EColi_output_$NAME\

bowtie2 -k 4 -p 6 -X2000 --local --mm -x $genomeIndexName -1 $INPUT_R1 -2 $INPUT_R2 --seed 23 2>$OUTPUT\.log | samtools view -Su /dev/stdin | samtools sort -n -o $OUTPUT\bowtie_EColi_output\.bam -

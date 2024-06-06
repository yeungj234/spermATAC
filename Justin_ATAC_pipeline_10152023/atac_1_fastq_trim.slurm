#! /bin/bash
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --verbose
#SBATCH -J cutadapt
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# example of sbatch command to run this script:
# sbatch -p risc --array=1-37 --export=inputfiles_R1=/lustre/fs4/risc_lab/scratch/jyeung/LS_PDvsDoxo_2/ATACseq/3022023_run_Deep_Justin_scripts/inputfiles_R1.txt,inputfiles_R2=/lustre/fs4/risc_lab/scratch/jyeung/LS_PDvsDoxo_2/ATACseq/3022023_run_Deep_Justin_scripts/inputfiles_R2.txt,outdir=/lustre/fs4/risc_lab/scratch/jyeung/LS_PDvsDoxo_2/ATACseq/3022023_run_Deep_Justin_scripts,samplelist=/lustre/fs4/risc_lab/scratch/jyeung/LS_PDvsDoxo_2/ATACseq/3022023_run_Deep_Justin_scripts/samplenames.txt /lustre/fs4/risc_lab/scratch/jyeung/scripts/pipelines/Justin_ATACpipeline/atac_1_fastq_trim.slurm

mkdir $outdir\/atac_1_fastq_trim/
INPUT=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID $samplelist | tail -n 1)
OUTPUT=$outdir\/atac_1_fastq_trim/$NAME\

source activate fastq2bam
cutadapt -a CTGTCTCTTATA -j 2 -e 0.1 -q 30 -m 5 -o $OUTPUT_\_trimmed.fastq $INPUT

fastqc $OUTPUT_\_trimmed.fastq


#! /bin/bash
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --verbose
#SBATCH -J cutadapt
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# example of sbatch command to run this script:
## sbatch -p risc --array=1-6 --export=inputfiles=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/Justin_ATAC_pipeline_10152023/atac_1_fastq_trim_interleaved/inputfiles.txt,outdir=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/Justin_ATAC_pipeline_10152023/atac_1_fastq_trim_interleaved,samplelist=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/Justin_ATAC_pipeline_10152023/atac_1_fastq_trim_interleaved/samples.txt /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/Justin_ATAC_pipeline_10152023/atac_1_fastq_trim_interleaved/atac_1_fastq_trim.slurm
mkdir $outdir\/atac_1_fastq_trim/
INPUT=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID $samplelist | tail -n 1)
OUTPUT=$outdir\/atac_1_fastq_trim/$NAME\

source activate fastq2bam
cutadapt --interleaved -a CTGTCTCTTATA -A CTGTCTCTTATA -j 2 -e 0.1 -q 30 -m 5 -o $OUTPUT_\_trimmed.fastq $INPUT

fastqc $OUTPUT_\_trimmed.fastq


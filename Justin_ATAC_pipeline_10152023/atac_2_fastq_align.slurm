#! /bin/bash
#SBATCH --verbose
#SBATCH -J bowtie2
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# example of sbatch command to run this script:
# sbatch -p risc --array=1-37 --export=outdir=/lustre/fs4/risc_lab/scratch/jyeung/LS_PDvsDoxo_2/ATACseq/3022023_run_Deep_Justin_scripts,genomeIndexName=/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/genome/Sequence/Bowtie2Index/genome,samplelist=/lustre/fs4/risc_lab/scratch/jyeung/LS_PDvsDoxo_2/ATACseq/3022023_run_Deep_Justin_scripts/samplenames.txt /lustre/fs4/risc_lab/scratch/jyeung/scripts/pipelines/Justin_ATACpipeline/atac_2_fastq_align.slurm

# make output directory for fastq align step
mkdir $outdir\/atac_2_fastq_align/

INPUT=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles_R1 | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID $samplelist | tail -n 1)
OUTPUT=$outdir\/atac_2_fastq_align/$NAME

source activate fastq2bam

bowtie2 -k 4 -p 6 -X2000 --local --mm -x $genomeIndexName -1 --interleaved $INPUT --seed 23 2>$OUTPUT\_$SLURM_JOB_ID\.log | samtools view -Su /dev/stdin | samtools sort -n -o $OUTPUT\_$SLURM_JOB_ID\.bam -

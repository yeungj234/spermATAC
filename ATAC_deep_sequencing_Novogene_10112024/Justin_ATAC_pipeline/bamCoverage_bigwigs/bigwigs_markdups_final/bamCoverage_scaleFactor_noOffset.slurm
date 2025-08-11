#! /bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --verbose
#SBATCH -J scaleFactor_bw
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

INPUT=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID $samplelist | tail -n 1)
SCALE=$(head -n $SLURM_ARRAY_TASK_ID $scalefactors | tail -n 1)
OUTPUT=${outdir}/${NAME}

# sbatch -p hpc_l40_b --array=1-6 --export=inputfiles=allreads_inputfiles.txt,outdir=/lustre/fs4/fnbk_lab/scratch/hkonishi/ATAC_deep_sequencing_Novogene_10112024/Justin_ATAC_pipeline/bamCoverage_bigwigs/bigwigs_markdups_final,samplelist=samplenames.txt,scalefactors=scaleFactor.txt bamCoverage_scaleFactor_noOffset.slurm

# keep entire fragment length, normalize by scaleFactor determined by DESeq2 and convert to bigwig for coverage plots. 

source activate ATACseq

# convert bam to bigwig using userspecified scale factors
bamCoverage -b $INPUT -o "${OUTPUT}_markdups_final_scaleFactornorm.bw" \
    --binSize 10 \
    --scaleFactor $SCALE \
    --extendReads \
    --numberOfProcessors 16

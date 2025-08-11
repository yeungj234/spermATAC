#!/bin/bash
#SBATCH --job-name=bamCoverage         # Job name
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8                     # Use 8 CPUs for parallelism
#SBATCH --mem=64G                             # 64 GB of memory
#SBATCH --time=24:00:00                       # Time limit (24 hours)
#SBATCH --output=%x-%j.out                    # Standard output
#SBATCH --error=%x-%j.err                     # Standard error
#SBATCH --verbose


INPUT=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID $samplelist | tail -n 1)
SCALE=$(head -n $SLURM_ARRAY_TASK_ID $scalefactors | tail -n 1)
OUTPUT=$outdir\/$NAME

# sbatch --array=1-10 --export=inputfiles=inputbams.txt,outdir=/lustre/fs4/risc_lab/scratch/jyeung/LS_PDvsDoxo_2/ATACseq/Tn5_offset_merged_RPGC_bigwigs/Tn5_offset_scaleFactornorm_bigwigs,samplelist=names.txt,scalefactors=scalefactor.txt bamCoverage_scaleFactor_noOffset.slurm
cd $outdir

# keep entire fragment length, normalize by scaleFactor determined by DESeq2 and convert to bigwig for coverage plots. 

source activate ATACseq
#samtools sort -o /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Justin_ATAC_pipeline/atac_2_fastq_align_05272025/$NAME\/$NAME\_sorted.bam $INPUT
#samtools index /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Justin_ATAC_pipeline/atac_2_fastq_align_05272025/$NAME\/$NAME\_sorted.bam
# convert bam to bigwig
bamCoverage -b /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Justin_ATAC_pipeline/atac_2_fastq_align_05272025/$NAME\/$NAME\_sorted.bam -o $OUTPUT\_scaleFactornorm.bw \
    --binSize 50 \
    --scaleFactor $SCALE \
    --extendReads \
    --numberOfProcessors 8

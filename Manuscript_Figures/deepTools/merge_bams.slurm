#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --verbose
#SBATCH -J mergebams
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# Submit with:
# sbatch -p hpc_bigmem_b /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/Manuscript_Figures/deepTools/merge_bams.slurm

########################################
# Activate conda environment
########################################

source $(conda info --base)/etc/profile.d/conda.sh
conda activate fastq2bam

########################################
# Define directories
########################################

BASE=/lustre/fs4/fnbk_lab/scratch/hkonishi/spermATAC/ATAC_deep_sequencing_Novogene_10112024/Justin_ATAC_pipeline/atac_4_bam_markduplicates

OUTDIR=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Justin_ATAC_pipeline/atac_4_bam_markduplicates/mergedbams

TMPDIR=$OUTDIR/tmp

mkdir -p $OUTDIR
mkdir -p $TMPDIR

########################################
# Merge Mitosis BAMs
########################################

samtools merge -@ 16 \
$OUTDIR/Mitosis_merged_markdups_final_multimap.bam \
$BASE/atacHAK002/atacHAK002_markdups_final.bam \
$BASE/atacHAK024/atacHAK024_markdups_final.bam \
$BASE/atacHAK025/atacHAK025_markdups_final.bam

########################################
# Sort Mitosis BAM
########################################

samtools sort -@ 8 -m 12G \
-T $TMPDIR/Mitosis \
-o $OUTDIR/Mitosis_merged_markdups_final_multimap.sorted.bam \
$OUTDIR/Mitosis_merged_markdups_final_multimap.bam

########################################
# Remove unsorted BAM
########################################

rm $OUTDIR/Mitosis_merged_markdups_final_multimap.bam

########################################
# Index sorted Mitosis BAM
########################################

samtools index -@ 8 \
$OUTDIR/Mitosis_merged_markdups_final_multimap.sorted.bam

########################################
# Merge Interphase BAMs
########################################

samtools merge -@ 16 \
$OUTDIR/Interphase_merged_markdups_final_multimap.bam \
$BASE/atacHAK001/atacHAK001_markdups_final.bam \
$BASE/atacHAK007/atacHAK007_markdups_final.bam \
$BASE/atacHAK023/atacHAK023_markdups_final.bam

########################################
# Sort Interphase BAM
########################################

samtools sort -@ 8 -m 12G \
-T $TMPDIR/Interphase \
-o $OUTDIR/Interphase_merged_markdups_final_multimap.sorted.bam \
$OUTDIR/Interphase_merged_markdups_final_multimap.bam

########################################
# Remove unsorted BAM
########################################

rm $OUTDIR/Interphase_merged_markdups_final_multimap.bam

########################################
# Index sorted Interphase BAM
########################################

samtools index -@ 8 \
$OUTDIR/Interphase_merged_markdups_final_multimap.sorted.bam
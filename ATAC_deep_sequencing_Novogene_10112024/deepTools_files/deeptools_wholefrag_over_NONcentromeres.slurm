#!/bin/bash
#SBATCH --job-name=deeptools_Coverage         # Job name
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8                     # Use 8 CPUs for parallelism                         # 64 GB of memory
#SBATCH --time=24:00:00                       # Time limit (24 hours)
#SBATCH --output=%x-%j.out                    # Standard output
#SBATCH --error=%x-%j.err                     # Standard error
#SBATCH --verbose

source activate encode-atac-seq-pipeline

BIGWIG_FILES="/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/scaleFactornorm_bigwigs/atacHAK001_merged_scaleFactornorm.bw
/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/scaleFactornorm_bigwigs/atacHAK007_merged_scaleFactornorm.bw
/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/scaleFactornorm_bigwigs/atacHAK023_merged_scaleFactornorm.bw
/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/scaleFactornorm_bigwigs/atacHAK002_merged_scaleFactornorm.bw
/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/scaleFactornorm_bigwigs/atacHAK024_merged_scaleFactornorm.bw
/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/scaleFactornorm_bigwigs/atacHAK025_merged_scaleFactornorm.bw"
samplenames="Interphase_1 Interphase_2 Interphase_3 Mitotic_1 Mitotic_2 Mitotic_3"
cd /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/deepTools_files
GENOMIC_REGIONS=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/peaks/nonCentromere_bins_based_on_GSE153058.bed
NAME="IvsM_ATAC_over_Non_Centromere_1kb_coverage"
computeMatrix scale-regions \
	-S $BIGWIG_FILES \
	-R $GENOMIC_REGIONS \
	-b 1000 -a 1000 -bs=10 -p=max \
	--regionBodyLength 1000 \
	--startLabel "start" \
	--endLabel "end" \
 	-out $NAME\_matrix.gz \
 	--numberOfProcessors "max" \
 	--skipZeros \
        --outFileSortedRegions $NAME\_sorted_regions.bed 

plotProfile -m $NAME\_matrix.gz \
    -out /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/deepTools_files/plotProfile/$NAME\_profile.pdf \
    --perGroup \
    --colors "#BADEFAFF" "#90CAF8FF" "#64B4F6FF" "#9FA7D9FF" "#5B6BBFFF" "#19227EFF" \
    --samplesLabel $samplenames \
	--startLabel "start" \
	--endLabel "end" \
    --plotTitle "$NAME" \
    --yMin 0 \
    --yMax 50 

plotHeatmap -m $NAME\_matrix.gz \
    -out  /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/deepTools_files/plotHeatmap/$NAME\_heatmap.pdf \
    --colorMap inferno \
    --samplesLabel $samplenames \
	--startLabel "start" \
	--endLabel "end" \
    --plotTitle "$NAME" \
    --yMin 0 \
    --yMax 50 
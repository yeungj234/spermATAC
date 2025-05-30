#!/bin/bash
#SBATCH --job-name=deeptools_Coverage         # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8                     # Use 8 CPUs for parallelism
#SBATCH --mem=64G                             # 64 GB of memory
#SBATCH --time=24:00:00                       # Time limit (24 hours)
#SBATCH --output=%x-%j.out                    # Standard output
#SBATCH --error=%x-%j.err                     # Standard error
#SBATCH --verbose

source activate encode-atac-seq-pipeline

mergedcuts="/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Tn5_offset_cutsCoverage/atacHAK001_scaleFactornorm_cutsCoverage.bw
/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Tn5_offset_cutsCoverage/atacHAK007_scaleFactornorm_cutsCoverage.bw
/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Tn5_offset_cutsCoverage/atacHAK023_scaleFactornorm_cutsCoverage.bw
/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Tn5_offset_cutsCoverage/atacHAK002_scaleFactornorm_cutsCoverage.bw
/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Tn5_offset_cutsCoverage/atacHAK024_scaleFactornorm_cutsCoverage.bw
/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Tn5_offset_cutsCoverage/atacHAK025_scaleFactornorm_cutsCoverage.bw"
samplenames="Interphase_1 Interphase_2 Interphase_3 Mitotic_1 Mitotic_2 Mitotic_3"
cd /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/deepTools_files
BIGWIG_FILES=$mergedcuts
GENOMIC_REGIONS=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/peaks/GSE153058_xla_v10.2_cen.bed
NAME="IvsM_ATAC_over_Centromere_5kb_coverage"
computeMatrix scale-regions \
	-S $BIGWIG_FILES \
	-R $GENOMIC_REGIONS \
	-b 5000 -a 5000 -bs=10 -p=max \
	--regionBodyLength 5000 \
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

plotHeatmap -m $NAME\_matrix.gz \
    -out  /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/deepTools_files/plotHeatmap/$NAME\_heatmap.pdf \
    --colorMap inferno \
    --samplesLabel $samplenames \
	--startLabel "start" \
	--endLabel "end" \
    --plotTitle "$NAME" \
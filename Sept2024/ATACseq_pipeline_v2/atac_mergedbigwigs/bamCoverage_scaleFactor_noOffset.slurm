#! /bin/bash
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --verbose
#SBATCH -J scaleFactor_bw
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# sbatch bamCoverage_scaleFactor_noOffset.slurm

source activate ATACseq

# blacklist_file=/rugpfs/fs0/risc_lab/scratch/jrendleman/ATAC_09262021/references/ENCFF940NTE_blacklist.bed
# convert bam to bigwig using userspecified scale factors
bamCoverage -b /lustre/fs4/risc_lab/store/jyeung/for_Hide/Sep2024/ATACseq_pipeline_v2/atac_mergedbams/atacHAK_Interphase_merged_sorted.bam -o /lustre/fs4/risc_lab/store/jyeung/for_Hide/Sep2024/ATACseq_pipeline_v2/atac_mergedbigwigs/atacHAK_Interphase_merged_scaleFactornorm.bw \
    --binSize 5 \
    --scaleFactor 0.9531843 \
    --extendReads \
    --numberOfProcessors 4
    
bamCoverage -b /lustre/fs4/risc_lab/store/jyeung/for_Hide/Sep2024/ATACseq_pipeline_v2/atac_mergedbams/atacHAK_Mitotic_merged_sorted.bam -o /lustre/fs4/risc_lab/store/jyeung/for_Hide/Sep2024/ATACseq_pipeline_v2/atac_mergedbigwigs/atacHAK_Mitotic_merged_scaleFactornorm.bw \
    --binSize 5 \
    --scaleFactor 0.7830034 \
    --extendReads \
    --numberOfProcessors 4

#! /bin/bash
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --verbose
#SBATCH -J cutsCoverageR
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

source /ru-auth/local/home/risc_soft/miniconda3/etc/profile.d/conda.sh
conda activate joanna

Rscript --vanilla /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Tn5_offset_cutsCoverage/cutsCoverage.R

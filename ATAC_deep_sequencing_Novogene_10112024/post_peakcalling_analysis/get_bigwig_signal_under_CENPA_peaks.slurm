#! /bin/bash
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --verbose
#SBATCH -J signal_in_peaks
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

source /ru-auth/local/home/risc_soft/miniconda3/etc/profile.d/conda.sh
conda activate joanna

Rscript --vanilla /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/post_peakcalling_analysis/get_bigwig_signal_under_CENPA_peaks.R
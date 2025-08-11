#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --verbose
#SBATCH -J picard.FLD
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# example of sbatch command to run this script:
# sbatch --array=1-9%3 --export=CENP_input=centromere_inputfiles.txt,NONCENP_input=noncentromere_inputfiles.txt,inputfiles=allreads_inputfiles.txt,outdir=$(pwd),samplelist=samplenames.txt atac_5_bam_FLD.slurm

INPUT_CENP=$(head -n $SLURM_ARRAY_TASK_ID $CENP_input | tail -n 1)
INPUT_NONCENP=$(head -n $SLURM_ARRAY_TASK_ID $NONCENP_input | tail -n 1)
INPUT_ALL=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID $samplelist | tail -n 1)
OUTPUT=$outdir\/$NAME

source activate fastq2bam

picard CollectInsertSizeMetrics INPUT=$INPUT_CENP OUTPUT=$OUTPUT\_centromeric_insertdata.txt H=$OUTPUT\_centromeric_FLD.pdf VERBOSITY=ERROR QUIET=TRUE USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE W=1000 
picard CollectInsertSizeMetrics INPUT=$INPUT_NONCENP OUTPUT=$OUTPUT\_noncentromeric_insertdata.txt H=$OUTPUT\_noncentromeric_FLD.pdf VERBOSITY=ERROR QUIET=TRUE USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE W=1000 
picard CollectInsertSizeMetrics INPUT=$INPUT_ALL OUTPUT=$OUTPUT\_allreads_insertdata.txt H=$OUTPUT\_allreads_FLD.pdf VERBOSITY=ERROR QUIET=TRUE USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE W=1000 
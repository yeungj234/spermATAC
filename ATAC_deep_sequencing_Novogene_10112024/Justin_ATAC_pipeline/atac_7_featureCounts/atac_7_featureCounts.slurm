#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --verbose
#SBATCH -J featureCounts
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# example of sbatch command to run this script:
# sbatch -p hpc_l40s_b --array=1-6 --export=inputfiles=allreads_inputfiles.txt,peaks=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Justin_ATAC_pipeline/atac_6_IDR_tagalign/nonredundant_masterpeaks_IDR.900000.0.1.gtf,outdir=$(pwd),samplelist=samplenames.txt atac_7_featureCounts.slurm

INPUT=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID $samplelist | tail -n 1)
# Define and create sample output directory
OUTPUT_DIR=${outdir}/${NAME}
mkdir -p "$OUTPUT_DIR"

source activate Hide_ATAC

# featureCounts -T 16 -p --countReadPairs -a $peaks -t peak -g peak_id -M --fraction -o "${OUTPUT_DIR}/${NAME}_counts.txt" $INPUT

# if an alignment overlaps more than 1 peak, assign it to the peak with the largest overlap to avoid counting a single alignment multiple times and avoid getting unassigned ambiguous reads. 
# featureCounts -T 16 -p --countReadPairs -a $peaks -t peak -g peak_id -M --largestOverlap --fraction -o "${OUTPUT_DIR}/${NAME}_largestOverlap_counts.txt" $INPUT

# still getting a few ambiguous reads assigned although less than before so try assigning reads to all overlaps instead of largestOverlap. 
featureCounts -T 16 -p --countReadPairs -a $peaks -t peak -g peak_id -M -O --fraction -o "${OUTPUT_DIR}/${NAME}_allOverlaps_counts.txt" $INPUT
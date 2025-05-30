#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --job-name=jf_kmer_query
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=jf_kmer_query-%j.out
#SBATCH --error=jf_kmer_query-%j.err
#SBATCH --verbose

# Usage: sbatch --array=1-6 --export=inputfiles_R1=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Centromeric_Repeat_Seq/inputfiles_R1.txt,inputfiles_R2=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Centromeric_Repeat_Seq/inputfiles_R2.txt,samplelist=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Centromeric_Repeat_Seq/samplenames.txt,outdir=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Centromeric_Repeat_Seq/jellyfish_query_output,MONOMER_KMERS=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Centromeric_Repeat_Seq/FCR_monomers/FCR_25mers_seqonly.txt jellyfish_kmer_query.slurm
READ1_GZ=$(head -n $SLURM_ARRAY_TASK_ID "$inputfiles_R1" | tail -n 1)
READ2_GZ=$(head -n $SLURM_ARRAY_TASK_ID "$inputfiles_R2" | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID "$samplelist" | tail -n 1)
OUTDIR="${outdir}/${NAME}"

source activate Hide_ATAC

# Output file paths
COUNT_DB_R1="${OUTDIR}_R1_atac_25mers.jf"
COUNT_DB_R2="${OUTDIR}_R2_atac_25mers.jf"
QUERY_OUT_R1="${OUTDIR}_R1_FCR_in_ATAC_counts.txt"
QUERY_OUT_R2="${OUTDIR}_R2_FCR_in_ATAC_counts.txt"

# Step 2: Count all 25-mers in FASTQ files
echo "Counting 25-mers in $READ1_GZ..."
jellyfish count -m 25 -s 30G -t 64 -C -o "$COUNT_DB_R1" <(zcat "$READ1_GZ")

echo "Counting 25-mers in $READ2_GZ..."
jellyfish count -m 25 -s 30G -t 64 -C -o "$COUNT_DB_R2" <(zcat "$READ2_GZ")

# Step 3: Query only FCR-derived 25-mers
echo "Querying FCR 25-mers against FASTQ counts..."
jellyfish query "$COUNT_DB_R1" < "$MONOMER_KMERS" > "$QUERY_OUT_R1"
jellyfish query "$COUNT_DB_R2" < "$MONOMER_KMERS" > "$QUERY_OUT_R2"

echo "Done. Results written to:"
echo "$QUERY_OUT_R1"
echo "$QUERY_OUT_R2"

#!/bin/bash
#SBATCH --ntasks 1                   
#SBATCH --cpus-per-task 2               
#SBATCH --time 4:00:00  
#SBATCH --mem-per-cpu 4G    
#SBATCH --verbose
#SBATCH -J fuzz_nuc
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# Required environment variables passed externally:
# $inputfiles_R1, $inputfiles_R2, $samplelist
# sbatch -p hpc_a10_a --array=1-6 --export=inputfiles_R1=inputfiles_R1.txt,inputfiles_R2=inputfiles_R2.txt,samplelist=samplenames.txt,outdir=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Centromeric_Repeat_Seq/fuzz_nuc_output/pmismatch_20,mismatch=20 fuzznuc_alignment.slurm

READ1_GZ=$(head -n $SLURM_ARRAY_TASK_ID "$inputfiles_R1" | tail -n 1)
READ2_GZ=$(head -n $SLURM_ARRAY_TASK_ID "$inputfiles_R2" | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID "$samplelist" | tail -n 1)

MONOMER_DIR="/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Centromeric_Repeat_Seq/FCR_monomers"

source activate Hide_ATAC

cd $outdir

echo "Decompressing FASTQ files..."
gunzip -c "$READ1_GZ" > $NAME\_tmp_R1.fastq
gunzip -c "$READ2_GZ" > $NAME\_tmp_R2.fastq

for monomer in "$MONOMER_DIR"/FCR_monomer_*.fasta; do
    base=$(basename "$monomer" .fasta)

    echo "Searching for $base in R1..."
    fuzznuc -sequence $NAME\_tmp_R1.fastq -pattern "$monomer" -pmismatch $mismatch -outfile "${outdir}/${NAME}_${base}_in_R1.txt"

    echo "Searching for $base in R2..."
    fuzznuc -sequence $NAME\_tmp_R2.fastq -pattern "$monomer" -pmismatch $mismatch -outfile "${outdir}/${NAME}_${base}_in_R2.txt"
done

echo "Cleaning up temporary files..."
rm $NAME\_tmp_R1.fastq $NAME\_tmp_R2.fastq


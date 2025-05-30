#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1                   
#SBATCH --cpus-per-task=16               # More CPUs for parallel execution
#SBATCH --time=24:00:00  
#SBATCH --mem=128G                       # Increase total memory for large FASTA input
#SBATCH --verbose
#SBATCH -J fuzz_nuc
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# Input variables are passed via --export
# Required: inputfiles_R1, inputfiles_R2, samplelist, MONOMER_KMERS, outdir, mismatch

READ1_GZ=$(head -n $SLURM_ARRAY_TASK_ID "$inputfiles_R1" | tail -n 1)
READ2_GZ=$(head -n $SLURM_ARRAY_TASK_ID "$inputfiles_R2" | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID "$samplelist" | tail -n 1)

source activate Hide_ATAC

# example sbatch command: 
# sbatch -p hpc_a10_a --array=1-6 --export=inputfiles_R1=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Centromeric_Repeat_Seq/inputfiles_R1.txt,inputfiles_R2=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Centromeric_Repeat_Seq/inputfiles_R2.txt,samplelist=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Centromeric_Repeat_Seq/samplenames.txt,outdir=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Centromeric_Repeat_Seq/fuzz_nuc_output/25mer,mismatch=2,MONOMER_KMERS=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Centromeric_Repeat_Seq/FCR_monomers/FCR_25mers_seqonly.fasta fuzznuc_kmer_alignment.slurm

cd "$outdir"

echo "[$(date)] Decompressing FASTQ files..."
seqtk seq -A "$READ1_GZ" > "${NAME}_tmp_R1.fasta" &
seqtk seq -A "$READ2_GZ" > "${NAME}_tmp_R2.fasta" &
wait

echo "[$(date)] Starting parallel fuzznuc searches..."

# Run R1 and R2 fuzznuc in parallel
fuzznuc -sequence "${NAME}_tmp_R1.fasta" -pattern "@$MONOMER_KMERS" -pmismatch "$mismatch" \
  -outfile "${outdir}/${NAME}_fuzznuc_output_in_R1.txt" &

fuzznuc -sequence "${NAME}_tmp_R2.fasta" -pattern "@$MONOMER_KMERS" -pmismatch "$mismatch" \
  -outfile "${outdir}/${NAME}_fuzznuc_output_in_R2.txt" &

wait
echo "[$(date)] fuzznuc search completed."

# Process R1
# split into chunk files to prevent memory allocation issues with grep. use grep to remove lines that start with # [0-9] with [0-9] being any number up to 6 digits
split -l 1000000 "${outdir}/${NAME}_fuzznuc_output_in_R1.txt" "${NAME}_R1_chunk_"
ls ${NAME}_R1_chunk_* | xargs -n 1 -P 4 -I{} bash -c '
  file={}
  grep -vE "^# [0-9]{1,6}" "$file" | grep -vE "^#=|^#-" > "'"${outdir}/${NAME}_"'${file}_R1_kmer_ids.txt"
'

# re-concatenate filtered chunk files: contains number of kmer hits per sequence and ids of kmers found in each sequence. 
cat "${outdir}/${NAME}_R1_chunk_"*_R1_kmer_ids.txt > "${outdir}/${NAME}_R1_kmer_ids_master.txt"
# filter for only Sequence ID from ATAC-seq files & number of hits corresponding to that sequence saved in separate file
cat "${outdir}/${NAME}_R1_kmer_ids_master.txt" |egrep 'Sequence|HitCount' | egrep "^#" > "${outdir}/${NAME}_R1_kmer_occurences.txt"
# filter and count the number of unique kmer occurrences in the ATAC-seq fasta file. 
cat "${outdir}/${NAME}_R1_kmer_ids_master.txt" |egrep -v "^#" |awk '{print $4}' |sed 's/^[^:]*://g' |sort |uniq -c |sort -nr >  "${outdir}/${NAME}_R1_kmer_ids_occurences.txt"

# Process R2
# split into chunk files to prevent memory allocation issues with grep. use grep to remove lines that start with # [0-9] with [0-9] being any number up to 6 digits

split -l 1000000 "${outdir}/${NAME}_fuzznuc_output_in_R2.txt" "${NAME}_R2_chunk_"
ls ${NAME}_R2_chunk_* | xargs -n 1 -P 4 -I{} bash -c '
  file={}
  grep -vE "^# [0-9]{1,6}" "$file" | grep -vE "^#=|^#-" > "'"${outdir}/${NAME}_"'${file}_R2_kmer_ids.txt"
'

# re-concatenate filtered chunk files: contains number of kmer hits per sequence and ids of kmers found in each sequence.
cat "${outdir}/${NAME}_R2_chunk_"*_R2_kmer_ids.txt > "${outdir}/${NAME}_R2_kmer_ids_master.txt"

# filter for only Sequence ID from ATAC-seq files & number of hits corresponding to that sequence saved in separate file
cat "${outdir}/${NAME}_R2_kmer_ids_master.txt" |egrep 'Sequence|HitCount' | egrep "^#" > "${outdir}/${NAME}_R2_kmer_occurences.txt"
# filter and count the number of unique kmer occurrences in the ATAC-seq fasta file. 
cat "${outdir}/${NAME}_R2_kmer_ids_master.txt" |egrep -v "^#" |awk '{print $4}' |sed 's/^[^:]*://g' |sort |uniq -c |sort -nr >  "${outdir}/${NAME}_R2_kmer_ids_occurences.txt"


echo "[$(date)] Cleaning up temporary files..."
rm "${NAME}_tmp_R1.fasta" "${NAME}_tmp_R2.fasta"
rm "${NAME}_R1_chunk_"* "${NAME}_R2_chunk_"*

echo "[$(date)] Done."

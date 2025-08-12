#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --verbose
#SBATCH -J merge

# example of sbatch command to run this script:
# sbatch --export=outdir=$(pwd),pval=0.05 bedtools_merge.slurm

cd $outdir
source activate fastq2bam

# make masterpeak set 
cat $(find |egrep '_pool.tn5_peaks.narrowPeak') > masterpeaks.narrowPeak

# Merge overlapping peaks 
bedtools merge -i masterpeaks.narrowPeak > nonredundant_masterpeaks.bed

# convert nonredundant master peak set to gtf for acceptable format for featureCounts to accept. 
awk 'BEGIN {OFS="\t"} 
{
  chrom=$1; start=$2+1; end=$3; name=chrom ":" start "-" end; score="."; strand="."; 
  print chrom, "narrowPeak", "peak", start, end, score, strand, ".", "peak_id \""name"\""
}' nonredundant_masterpeaks.bed > nonredundant_masterpeaks.gtf

# make masterpeak set after IDR filtering
zcat $(find $outdir | egrep "_REP1_VS_REP2_VS_REP3\.IDR.${pval}.narrowPeak\.gz$") > masterpeaks_IDR.${pval}.narrowPeak
sort -k1,1 -k2,2n masterpeaks_IDR.${pval}.narrowPeak > masterpeaks_IDR\.${pval}\.sorted\.narrowPeak

# merge overlapping IDR peaks
bedtools merge -i masterpeaks_IDR\.${pval}\.sorted\.narrowPeak > nonredundant_masterpeaks_IDR\.${pval}\.bed

# convert nonredundant master peak set to gtf for acceptable format for featureCounts to accept. 
awk 'BEGIN {OFS="\t"} 
{
  chrom=$1; start=$2+1; end=$3; name=chrom ":" start "-" end; score="."; strand="."; 
  print chrom, "narrowPeak", "peak", start, end, score, strand, ".", "peak_id \""name"\""
}' nonredundant_masterpeaks_IDR\.${pval}\.bed > nonredundant_masterpeaks_IDR\.${pval}\.gtf

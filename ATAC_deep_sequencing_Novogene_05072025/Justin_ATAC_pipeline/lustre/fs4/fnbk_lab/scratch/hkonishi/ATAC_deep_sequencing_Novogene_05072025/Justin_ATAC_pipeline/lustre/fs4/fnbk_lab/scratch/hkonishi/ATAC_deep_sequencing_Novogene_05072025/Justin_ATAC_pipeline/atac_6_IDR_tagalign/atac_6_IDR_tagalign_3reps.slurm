#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=384G
#SBATCH --time=24:00:00
#SBATCH --verbose
#SBATCH -J tagAlign
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# example of sbatch command to run this script:
# sbatch -p hpc_a10_a --array=1-3 --export=inputfiles_Rep1=inputfiles_Rep1.txt,inputfiles_Rep2=inputfiles_Rep2.txt,inputfiles_Rep3=inputfiles_Rep3.txt,outdir=$(pwd),samplelist=condition_names.txt atac_6_IDR_tagalign_3reps.slurm

INPUT_Rep1=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles_Rep1 | tail -n 1)
INPUT_Rep2=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles_Rep2 | tail -n 1)
INPUT_Rep3=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles_Rep3 | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID $samplelist | tail -n 1)

# Define and create sample output directory
OUTPUT_DIR=${outdir}/${NAME}
mkdir -p "$OUTPUT_DIR"

source activate fastq2bam

#Rep 1
#samtools sort -@ 24 -n $INPUT_Rep1 | bedtools bamtobed -bedpe -mate1 -i - | gzip -nc > ${OUTPUT_DIR}/${NAME}_Rep1.bedpe.gz
#zcat ${OUTPUT_DIR}/${NAME}_Rep1.bedpe.gz | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | \gzip -nc > ${OUTPUT_DIR}/${NAME}_Rep1.tagAlign.gz

#Rep 2
#samtools sort -@ 24 -n $INPUT_Rep2 | bedtools bamtobed -bedpe -mate1 -i - | gzip -nc > ${OUTPUT_DIR}/${NAME}_Rep2.bedpe.gz
#zcat ${OUTPUT_DIR}/${NAME}_Rep2.bedpe.gz | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | \gzip -nc > ${OUTPUT_DIR}/${NAME}_Rep2.tagAlign.gz

#Rep 3
#samtools sort -@ 24 -n $INPUT_Rep3 | bedtools bamtobed -bedpe -mate1 -i - | gzip -nc > ${OUTPUT_DIR}/${NAME}_Rep3.bedpe.gz
#zcat ${OUTPUT_DIR}/${NAME}_Rep3.bedpe.gz | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | \gzip -nc > ${OUTPUT_DIR}/${NAME}_Rep3.tagAlign.gz

#zcat ${OUTPUT_DIR}/${NAME}_Rep1.tagAlign.gz ${OUTPUT_DIR}/${NAME}_Rep2.tagAlign.gz ${OUTPUT_DIR}/${NAME}_Rep3.tagAlign.gz | gzip -nc > ${OUTPUT_DIR}/${NAME}_pool.tagAlign.gz

#Tn5 shifting of TagAlign files to account for Tn5 insertion by +4bp on + strand and -5bp on - strand
#zcat ${OUTPUT_DIR}/${NAME}_Rep1.tagAlign.gz | awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6=="+") {$2=$2+4} else if ($6=="-") {$3=$3-5} print $0}' | gzip -nc > ${OUTPUT_DIR}/${NAME}_Rep1.tn5.tagAlign.gz
#zcat ${OUTPUT_DIR}/${NAME}_Rep2.tagAlign.gz | awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6=="+") {$2=$2+4} else if ($6=="-") {$3=$3-5} print $0}' | gzip -nc > ${OUTPUT_DIR}/${NAME}_Rep2.tn5.tagAlign.gz
#zcat ${OUTPUT_DIR}/${NAME}_Rep3.tagAlign.gz | awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6=="+") {$2=$2+4} else if ($6=="-") {$3=$3-5} print $0}' | gzip -nc > ${OUTPUT_DIR}/${NAME}_Rep3.tn5.tagAlign.gz
#zcat ${OUTPUT_DIR}/${NAME}_pool.tagAlign.gz | awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6=="+") {$2=$2+4} else if ($6=="-") {$3=$3-5} print $0}' | gzip -nc > ${OUTPUT_DIR}/${NAME}_pool.tn5.tagAlign.gz

conda deactivate
source activate macs2_python2

##Call peaks replicates
#macs2 callpeak -t ${OUTPUT_DIR}/${NAME}_Rep1.tn5.tagAlign.gz -f BED -n ${OUTPUT_DIR}/${NAME}_Rep1.tn5 -g 2.7e9 -p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits
#macs2 callpeak -t ${OUTPUT_DIR}/${NAME}_Rep2.tn5.tagAlign.gz -f BED -n ${OUTPUT_DIR}/${NAME}_Rep2.tn5 -g 2.7e9 -p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits
#macs2 callpeak -t ${OUTPUT_DIR}/${NAME}_Rep3.tn5.tagAlign.gz -f BED -n ${OUTPUT_DIR}/${NAME}_Rep3.tn5 -g 2.7e9 -p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits
##Call peaks pooled sample
macs2 callpeak -t ${OUTPUT_DIR}/${NAME}_pool.tn5.tagAlign.gz -f BED -n ${OUTPUT_DIR}/${NAME}_pool.tn5 -g 2.7e9 -p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits

# typically top 300,000 highest signal peaks but we will increase it to 900,000 peaks to keep more peaks for this specific kind of analysis. 
sort -k 8gr,8gr ${OUTPUT_DIR}/${NAME}_Rep1.tn5_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | head -n 900000 | gzip -nc > ${OUTPUT_DIR}/${NAME}_Rep1.tn5.narrowPeak.gz
sort -k 8gr,8gr ${OUTPUT_DIR}/${NAME}_Rep2.tn5_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | head -n 900000 | gzip -nc > ${OUTPUT_DIR}/${NAME}_Rep2.tn5.narrowPeak.gz
sort -k 8gr,8gr ${OUTPUT_DIR}/${NAME}_Rep3.tn5_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | head -n 900000 | gzip -nc > ${OUTPUT_DIR}/${NAME}_Rep3.tn5.narrowPeak.gz
sort -k 8gr,8gr ${OUTPUT_DIR}/${NAME}_pool.tn5_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | head -n 900000 | gzip -nc > ${OUTPUT_DIR}/${NAME}_pool.tn5.narrowPeak.gz

#apply IDR

conda deactivate
source activate encode-atac-seq-pipeline
idr --samples ${OUTPUT_DIR}/${NAME}_Rep1.tn5.narrowPeak.gz ${OUTPUT_DIR}/${NAME}_Rep2.tn5.narrowPeak.gz \
--peak-list ${OUTPUT_DIR}/${NAME}_pool.tn5.narrowPeak.gz --input-file-type narrowPeak \
--output-file ${OUTPUT_DIR}/${NAME}_IDR_12.tn5.narrowPeak --rank p.value --soft-idr-threshold 0.05 \
--plot --use-best-multisummit-IDR

idr --samples ${OUTPUT_DIR}/${NAME}_Rep2.tn5.narrowPeak.gz ${OUTPUT_DIR}/${NAME}_Rep3.tn5.narrowPeak.gz \
--peak-list ${OUTPUT_DIR}/${NAME}_pool.tn5.narrowPeak.gz --input-file-type narrowPeak \
--output-file ${OUTPUT_DIR}/${NAME}_IDR_23.tn5.narrowPeak --rank p.value --soft-idr-threshold 0.05 \
--plot --use-best-multisummit-IDR

idr --samples  ${OUTPUT_DIR}/${NAME}_IDR_12.tn5.narrowPeak  ${OUTPUT_DIR}/${NAME}_IDR_23.tn5.narrowPeak \
--peak-list ${OUTPUT_DIR}/${NAME}_pool.tn5.narrowPeak.gz --input-file-type narrowPeak \
--output-file ${OUTPUT_DIR}/${NAME}_IDR_final.tn5.narrowPeak --rank p.value --soft-idr-threshold 0.05 \
--plot --use-best-multisummit-IDR

IDR_THRESH_TRANSFORMED=$(awk -v p=0.05 'BEGIN{print -log(p)/log(10)}')

awk 'BEGIN{OFS="\t"} $12>='"${IDR_THRESH_TRANSFORMED}"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ${OUTPUT_DIR}/${NAME}_IDR_final.tn5.narrowPeak | \
sort | uniq | sort -k7n,7n | gzip -nc > ${OUTPUT_DIR}/${NAME}_REP1_VS_REP2_VS_REP3.IDR0.05.narrowPeak.gz

# skip blacklist filtering step because we don't have that for Xenopus + want to keep reads mapping to repetitive regions

# generate redundant master peak set. 
cat $(find -L $PWD -name "*_REP1_VS_REP2_VS_REP3.IDR0.05.narrowPeak.gz") > $outdir/master_IDR_900000.Peaks.narrowPeak.gz


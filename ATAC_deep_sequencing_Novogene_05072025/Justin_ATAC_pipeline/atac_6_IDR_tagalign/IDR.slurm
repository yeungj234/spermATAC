#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --verbose
#SBATCH -J IDR
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# example of sbatch command to run this script:
# sbatch -p hpc_l40_a --array=1-3 --export=outdir=$(pwd),pval=0.1,samplelist=condition_names.txt IDR.slurm

NAME=$(head -n $SLURM_ARRAY_TASK_ID $samplelist | tail -n 1)
OUTPUT_DIR=${outdir}/${NAME}

#source activate macs2_python2

# typically top 300,000 highest signal peaks but I am testing out a higher cut-off value of 600,000 to keep more peaks. 
#sort -k 8gr,8gr ${OUTPUT_DIR}/${NAME}_Rep1.tn5_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | head -n $peaknum | gzip -nc > ${OUTPUT_DIR}/${NAME}_Rep1.tn5.narrowPeak.gz
#sort -k 8gr,8gr ${OUTPUT_DIR}/${NAME}_Rep2.tn5_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | head -n $peaknum | gzip -nc > ${OUTPUT_DIR}/${NAME}_Rep2.tn5.narrowPeak.gz
#sort -k 8gr,8gr ${OUTPUT_DIR}/${NAME}_Rep3.tn5_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | head -n $peaknum | gzip -nc > ${OUTPUT_DIR}/${NAME}_Rep3.tn5.narrowPeak.gz
#sort -k 8gr,8gr ${OUTPUT_DIR}/${NAME}_pool.tn5_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | head -n $peaknum | gzip -nc > ${OUTPUT_DIR}/${NAME}_pool.tn5.narrowPeak.gz

#apply IDR

#conda deactivate
source activate encode-atac-seq-pipeline
idr --samples ${OUTPUT_DIR}/${NAME}_Rep1.tn5.narrowPeak.gz ${OUTPUT_DIR}/${NAME}_Rep2.tn5.narrowPeak.gz \
--peak-list ${OUTPUT_DIR}/${NAME}_pool.tn5.narrowPeak.gz --input-file-type narrowPeak \
--output-file ${OUTPUT_DIR}/${NAME}_IDR_12.${pval}.tn5.narrowPeak --rank p.value --soft-idr-threshold $pval \
--plot --use-best-multisummit-IDR

idr --samples ${OUTPUT_DIR}/${NAME}_Rep2.tn5.narrowPeak.gz ${OUTPUT_DIR}/${NAME}_Rep3.tn5.narrowPeak.gz \
--peak-list ${OUTPUT_DIR}/${NAME}_pool.tn5.narrowPeak.gz --input-file-type narrowPeak \
--output-file ${OUTPUT_DIR}/${NAME}_IDR_23.${pval}.tn5.narrowPeak --rank p.value --soft-idr-threshold $pval \
--plot --use-best-multisummit-IDR

idr --samples  ${OUTPUT_DIR}/${NAME}_IDR_12.${pval}.tn5.narrowPeak  ${OUTPUT_DIR}/${NAME}_IDR_23.${pval}.tn5.narrowPeak \
--peak-list ${OUTPUT_DIR}/${NAME}_pool.tn5.narrowPeak.gz --input-file-type narrowPeak \
--output-file ${OUTPUT_DIR}/${NAME}_IDR_final.${pval}.tn5.narrowPeak --rank p.value --soft-idr-threshold $pval \
--plot --use-best-multisummit-IDR

IDR_THRESH_TRANSFORMED=$(awk -v p=$pval 'BEGIN{print -log(p)/log(10)}')

awk 'BEGIN{OFS="\t"} $12>='"${IDR_THRESH_TRANSFORMED}"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ${OUTPUT_DIR}/${NAME}_IDR_final.${pval}.tn5.narrowPeak | \
sort | uniq | sort -k7n,7n | gzip -nc > ${OUTPUT_DIR}/${NAME}_REP1_VS_REP2_VS_REP3.IDR.${pval}.narrowPeak.gz

# skip blacklist filtering step because we don't have that for Xenopus + want to keep reads mapping to repetitive regions

# generate redundant master peak set. 
zcat $(find -L $PWD -name "*_REP1_VS_REP2_VS_REP3.IDR.${pval}.narrowPeak.gz") > $outdir/master_IDR.${pval}_Peaks.narrowPeak.gz



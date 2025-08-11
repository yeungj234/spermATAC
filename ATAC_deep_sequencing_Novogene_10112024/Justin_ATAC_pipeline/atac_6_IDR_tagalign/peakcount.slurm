#! /bin/bash
echo "Filename,LineCount" > peak_counts_masterpeaks.csv
ls /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Justin_ATAC_pipeline/atac_6_IDR_tagalign/ -type f \( -name "*.narrowPeak" -o -name "*.bed" \) | while read file; do
    count=$(wc -l < "$file")
    filename=$(basename "$file")
    echo "$filename,$count" >> peak_counts_masterpeaks.csv
done


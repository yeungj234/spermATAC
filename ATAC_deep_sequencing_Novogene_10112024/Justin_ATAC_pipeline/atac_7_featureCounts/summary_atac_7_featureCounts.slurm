#!/bin/bash

# Directory with flagstat files
featureCounts_dir=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Justin_ATAC_pipeline/atac_7_featureCounts

# Output CSV file
OUTPUT_CSV="${featureCounts_dir}/featureCounts_ReadsinPeaks.csv"

# files
cd $featureCounts_dir
flagstat_files=$(find $pwd -name '*_counts.txt.summary' |sort -n)
# Write header
echo "Sample,ReadsinPeaks,NonOverlappingReads,AmbiguousReads,FRIP" > "$OUTPUT_CSV"
# Loop through each flagstat file
for file in $flagstat_files; do
    sample=$(basename $file _counts.txt.summary)

    # Extract stats using awk
    read readsinpeaks _ < <(awk 'NR==2{print $2}' "$file")
    read readsNotinPeaks _ < <(awk 'NR==13{print $2}' "$file")
    read ambiguity _ < <(awk 'NR==15{print $2}' "$file")
    FRIP=$(awk -v a="$readsinpeaks" -v b="$readsNotinPeaks" -v c="$ambiguity" 'BEGIN {print a / (a + b + c)}')

    echo "$sample,$readsinpeaks,$readsNotinPeaks,$ambiguity,$FRIP" >> "$OUTPUT_CSV"
done

echo "Summary written to $OUTPUT_CSV"


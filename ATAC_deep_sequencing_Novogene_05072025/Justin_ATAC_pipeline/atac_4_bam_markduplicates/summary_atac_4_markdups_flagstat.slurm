#!/bin/bash

# Directory with flagstat files
FLAGSTAT_DIR="/lustre/fs4/fnbk_lab/scratch/hkonishi/spermATAC/ATAC_deep_sequencing_Novogene_05072025/Justin_ATAC_pipeline/atac_4_bam_markduplicates"

# Output CSV file
OUTPUT_CSV="${FLAGSTAT_DIR}/flagstat_atac_4_filter_summary.csv"

# files
cd $FLAGSTAT_DIR
flagstat_files=$(find $pwd -name '*_filter.flagstat' |sort -n)
# Write header
echo "Sample,TotalReads,SecondaryReads,SupplementaryReads,Duplicates,Mapped,Paired,Read1,Read2,ProperlyPaired,WithMateMapped,Singletons,MateDiffChr,MateDiffChrMapQover5" > "$OUTPUT_CSV"
# Loop through each flagstat file
for file in $flagstat_files; do
    sample=$(basename $file _filter.flagstat |sed 's/_[0-9]\{7\}$//')

    # Extract stats using awk
    read total _ < <(awk 'NR==1{print $1}' "$file")
    read secondary _ < <(awk 'NR==2{print $1}' "$file")
    read suppl _ < <(awk 'NR==3{print $1}' "$file")
    read dupes _ < <(awk 'NR==4{print $1}' "$file")
    read mapped _ < <(awk 'NR==5{print $1}' "$file")
    read paired _ < <(awk 'NR==6{print $1}' "$file")
    read read1 _ < <(awk 'NR==7{print $1}' "$file")
    read read2 _ < <(awk 'NR==8{print $1}' "$file")
    read properpair _ < <(awk 'NR==9{print $1}' "$file")
    read withMate _ < <(awk 'NR==10{print $1}' "$file")
    read singletons _ < <(awk 'NR==11{print $1}' "$file")
    read mateDiff _ < <(awk 'NR==12{print $1}' "$file")
    read mateDiffQ _ < <(awk 'NR==13{print $1}' "$file")

    echo "$sample,$total,$secondary,$suppl,$dupes,$mapped,$paired,$read1,$read2,$properpair,$withMate,$singletons,$mateDiff,$mateDiffQ" >> "$OUTPUT_CSV"
done

echo "Summary written to $OUTPUT_CSV"

# Output CSV file
OUTPUT_CSV="${FLAGSTAT_DIR}/flagstat_atac_4_fixmate_summary.csv"

# files
cd $FLAGSTAT_DIR
flagstat_files=$(find $pwd -name '*_fixmate.flagstat' |sort -n)
# Write header
echo "Sample,TotalReads,SecondaryReads,SupplementaryReads,Duplicates,Mapped,Paired,Read1,Read2,ProperlyPaired,WithMateMapped,Singletons,MateDiffChr,MateDiffChrMapQover5" > "$OUTPUT_CSV"
# Loop through each flagstat file
for file in $flagstat_files; do
    sample=$(basename $file _fixmate.flagstat |sed 's/_[0-9]\{7\}$//')

    # Extract stats using awk
    read total _ < <(awk 'NR==1{print $1}' "$file")
    read secondary _ < <(awk 'NR==2{print $1}' "$file")
    read suppl _ < <(awk 'NR==3{print $1}' "$file")
    read dupes _ < <(awk 'NR==4{print $1}' "$file")
    read mapped _ < <(awk 'NR==5{print $1}' "$file")
    read paired _ < <(awk 'NR==6{print $1}' "$file")
    read read1 _ < <(awk 'NR==7{print $1}' "$file")
    read read2 _ < <(awk 'NR==8{print $1}' "$file")
    read properpair _ < <(awk 'NR==9{print $1}' "$file")
    read withMate _ < <(awk 'NR==10{print $1}' "$file")
    read singletons _ < <(awk 'NR==11{print $1}' "$file")
    read mateDiff _ < <(awk 'NR==12{print $1}' "$file")
    read mateDiffQ _ < <(awk 'NR==13{print $1}' "$file")

    echo "$sample,$total,$secondary,$suppl,$dupes,$mapped,$paired,$read1,$read2,$properpair,$withMate,$singletons,$mateDiff,$mateDiffQ" >> "$OUTPUT_CSV"
done

echo "Summary written to $OUTPUT_CSV"

# Output CSV file
OUTPUT_CSV="${FLAGSTAT_DIR}/flagstat_atac_4_markdups_final_summary.csv"

# files
cd $FLAGSTAT_DIR
flagstat_files=$(find $pwd -name '*_markdups_final.flagstat' |sort -n)
# Write header
echo "Sample,TotalReads,SecondaryReads,SupplementaryReads,Duplicates,Mapped,Paired,Read1,Read2,ProperlyPaired,WithMateMapped,Singletons,MateDiffChr,MateDiffChrMapQover5" > "$OUTPUT_CSV"
# Loop through each flagstat file
for file in $flagstat_files; do
    sample=$(basename $file _markdups_final.flagstat |sed 's/_[0-9]\{7\}$//')

    # Extract stats using awk
    read total _ < <(awk 'NR==1{print $1}' "$file")
    read secondary _ < <(awk 'NR==2{print $1}' "$file")
    read suppl _ < <(awk 'NR==3{print $1}' "$file")
    read dupes _ < <(awk 'NR==4{print $1}' "$file")
    read mapped _ < <(awk 'NR==5{print $1}' "$file")
    read paired _ < <(awk 'NR==6{print $1}' "$file")
    read read1 _ < <(awk 'NR==7{print $1}' "$file")
    read read2 _ < <(awk 'NR==8{print $1}' "$file")
    read properpair _ < <(awk 'NR==9{print $1}' "$file")
    read withMate _ < <(awk 'NR==10{print $1}' "$file")
    read singletons _ < <(awk 'NR==11{print $1}' "$file")
    read mateDiff _ < <(awk 'NR==12{print $1}' "$file")
    read mateDiffQ _ < <(awk 'NR==13{print $1}' "$file")

    echo "$sample,$total,$secondary,$suppl,$dupes,$mapped,$paired,$read1,$read2,$properpair,$withMate,$singletons,$mateDiff,$mateDiffQ" >> "$OUTPUT_CSV"
done

echo "Summary written to $OUTPUT_CSV"
#!/bin/bash

# Directory with flagstat files
FLAGSTAT_DIR="/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_05072025/Justin_ATAC_pipeline/atac_2_fastq_align"

# Output CSV file
OUTPUT_CSV="${FLAGSTAT_DIR}/flagstat_atac_2_fastq_align_summary.csv"

# files
cd $FLAGSTAT_DIR
flagstat_files=$(find $pwd -name '*flagstat.txt' |sort -n)
# Write header
echo "Sample,TotalReads,SecondaryReads,SupplementaryReads,Duplicates,Mapped,Paired,Read1,Read2,ProperlyPaired,WithMateMapped,Singletons,MateDiffChr,MateDiffChrMapQover5" > "$OUTPUT_CSV"
# Loop through each flagstat file
for file in $flagstat_files; do
    sample=$(basename $file _flagstat.txt |sed 's/_[0-9]\{7\}$//')

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


# log files
log_files=$(find $pwd -name '*.log' |sort -n)
# Output CSV file
OUTPUT_LOG="${FLAGSTAT_DIR}/alignmentStats_atac_2_fastq_align_summary.csv"
# Write header
echo "Sample,TotalReads,Paired,concordance_0,concordance_1,concordance_over1,Percent_overallAlignment" > "$OUTPUT_LOG"
# Loop through each flagstat file
for file in $log_files; do
    sample=$(basename $file .log |sed 's/_[0-9]\{7\}$//')

    # Extract stats using awk
    read total _ < <(awk 'NR==1{print $1}' "$file")
    read paired _ < <(awk 'NR==2{print $1}' "$file")
    read concord0 _ < <(awk 'NR==3{print $1}' "$file")
    read concord1 _ < <(awk 'NR==4{print $1}' "$file")
    read concordover1 _ < <(awk 'NR==5{print $1}' "$file")
    read alignrate < <(awk 'NR==15{print $1}' "$file")

    echo "$sample,$total,$paired,$concord0,$concord1,$concordover1,$alignrate" >> "$OUTPUT_LOG"
done


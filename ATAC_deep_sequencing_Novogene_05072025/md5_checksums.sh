#!/bin/bash

# Output CSV file name
output_file="checksums.csv"
files=$(find /rugpfs/fs0/cem/store/hkonishi/Data_Novogene/May2025/atacSeq/01.RawData/ -name '*fq.gz' |sort)
# Print the header to the CSV file
echo "File Name,MD5 Checksum" > "$output_file"

# Loop through all files in the current directory
for file in $files; do
  # Check if the item is a file (not a directory)
  if [[ -f "$file" ]]; then
    # Compute the MD5 checksum
    checksum=$(md5sum "$file" | awk '{ print $1 }')
    # Append the file name and checksum to the CSV file
    echo "\"$file\",\"$checksum\"" >> "$output_file"
  fi
done

# Print completion message
echo "MD5 checksums have been saved to $output_file"

#!/bin/bash
# Loop over all files ending with _counts.tsv in the current directory
for f in *_counts.tsv; do
    # Generate an output filename by replacing .tsv with _summary.tsv
    out="${f%.tsv}_summary.tsv"

    # Begin command block:
    (
        # Print header line with column names
        echo -e "kmer_frequencies\tNumber_of_reads"

        # Process the file:
        cut -f2 "$f" |            # Extract the second column (k-mer frequency)
        sort |                    # Sort the values
        uniq -c |                 # Count occurrences of each unique value
        sort -nr |                # Sort numerically in reverse (so highest count comes first)
        sort -n -k2 |             # Then sort by k-mer frequency value (column 2)
        awk '{print $2 "\t" $1}'  # Swap columns: frequency first, then count
    ) > "$out"                    # Write the entire output to the new summary file
done


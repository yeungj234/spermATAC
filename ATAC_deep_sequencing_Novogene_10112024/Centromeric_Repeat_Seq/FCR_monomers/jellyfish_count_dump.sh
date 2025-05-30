#!/bin/bash
conda activate Hide_ATAC
# count 25mers
jellyfish count -m 25 -s 10M -t 4 -C -o fcr_25mers.jf GSE153058_FCR_monomers.fa
# dump sequences of 25mers (contains sequences and number of occurrences)
jellyfish dump -c fcr_25mers.jf > fcr_25mers.txt
# keep only sequences
cat fcr_25mers.txt |sed 's/ [0-9]//g' > FCR_25mers_seqonly.txt
# convert to fasta format. 
awk '{print ">" NR "\n" $0}' FCR_25mers_seqonly.txt > FCR_25mers_seqonly.fasta


# count 10-mers
jellyfish count -m 10 -s 10M -t 4 -C -o fcr_10mers.jf GSE153058_FCR_monomers.fa
# dump sequences of 10mers (contains sequences and number of occurrences)
jellyfish dump -c fcr_10mers.jf > fcr_10mers.txt
# keep only sequences
cat fcr_10mers.txt |sed 's/ [0-9]//g' > FCR_10mers_seqonly.txt
# convert to fasta format. 
awk '{print ">" NR "\n" $0}' FCR_10mers_seqonly.txt > FCR_10mers_seqonly.fasta

# count 15-mers
jellyfish count -m 15 -s 10M -t 4 -C -o fcr_15mers.jf GSE153058_FCR_monomers.fa
# dump sequences of 15mers (contains sequences and number of occurrences)
jellyfish dump -c fcr_15mers.jf > fcr_15mers.txt
# keep only sequences
cat fcr_15mers.txt |sed 's/ [0-9]//g' > FCR_15mers_seqonly.txt
# convert to fasta format. 
awk '{print ">" NR "\n" $0}' FCR_15mers_seqonly.txt > FCR_15mers_seqonly.fasta
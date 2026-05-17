import gzip  # Module to read compressed .gz files
import argparse  # Module to parse command-line arguments
import ahocorasick  # For building and using Aho-Corasick trie (efficient multi-pattern matching)
from Bio import SeqIO  # From Biopython, for parsing FASTA and FASTQ files

# Function: Load k-mers from a FASTA file (supports .gz compressed files)
def load_kmers_from_fasta(fasta_file):
    kmers = {}  # Dictionary to store {kmer_sequence: kmer_id}
    
    # Choose open function based on whether the file is gzipped
    open_func = gzip.open if fasta_file.endswith('.gz') else open
    
    with open_func(fasta_file, "rt") as handle:  # Open the file in text mode
        for record in SeqIO.parse(handle, "fasta"):  # Iterate through FASTA records
            kmers[str(record.seq)] = record.id  # Store kmer sequence as key, FASTA ID as value
    return kmers  # Return dictionary of all kmers

# Function: Build the Aho-Corasick automaton from the list of kmers
def build_automaton(kmers):
    A = ahocorasick.Automaton()  # Initialize an empty Aho-Corasick automaton
    
    # Add each kmer into the automaton
    for kmer_seq, kmer_id in kmers.items():
        A.add_word(kmer_seq, (kmer_id, kmer_seq))  # Store a tuple: (ID, sequence)
    
    A.make_automaton()  # Finalize the automaton (important!)
    return A  # Return the fully constructed automaton

# Function: Match kmers to FASTQ reads and write two output files
def count_kmers_in_fastq(fastq_file, automaton, output_count_file=None, output_detail_file=None):
    # Choose appropriate file opening method
    open_func = gzip.open if fastq_file.endswith(".gz") else open

    with open_func(fastq_file, "rt") as handle, \
         open(output_count_file, "w") as count_out, \
         open(output_detail_file, "w") as detail_out:

        # Write column headers for both outputs
        count_out.write("read_id\tnum_kmers_matched\n")
        detail_out.write("read_id\tmatched_kmer_ids\n")

        # Iterate over each FASTQ record
        for record in SeqIO.parse(handle, "fastq"):
            read_id = record.id  # Get read identifier (e.g. "@SEQ_ID")
            read_seq = str(record.seq)  # Extract read sequence as string
            matched_ids = set()  # Store matched kmer IDs as a set (unique values only)

            # Search for all matching kmers in the read
            for end_index, (kmer_id, kmer_seq) in automaton.iter(read_seq):
                matched_ids.add(kmer_id)  # Store only the kmer ID (not full sequence)

            # Count and output total matches
            num_matches = len(matched_ids)
            count_out.write(f"{read_id}\t{num_matches}\n")
            detail_out.write(f"{read_id}\t{','.join(matched_ids)}\n")  # Output list of matched kmer IDs

# This runs if the script is executed as a standalone program
if __name__ == "__main__":
    # Define command-line arguments
    parser = argparse.ArgumentParser(
        description="Count and identify k-mer matches per read using Aho-Corasick."
    )
    
    # Required: K-mer FASTA file
    parser.add_argument(
        "-k", "--kmers", required=True, help="FASTA file with kmers (supports .gz)"
    )

    # Required: Input FASTQ file (supports gzipped)
    parser.add_argument(
        "-f", "--fastq", required=True, help="FASTQ file with reads (supports .gz)"
    )

    # Output file for counts per read
    parser.add_argument(
        "-o", "--output_count", help="Output file for match counts per read"
    )

    # Output file for IDs of matched kmers per read
    parser.add_argument(
        "-d", "--output_detail", help="Output file listing matched kmers per read"
    )

    # Parse all arguments
    args = parser.parse_args()

    # Load kmers and build the Aho-Corasick automaton
    kmers = load_kmers_from_fasta(args.kmers)
    automaton = build_automaton(kmers)

    # Run kmer search on FASTQ file
    count_kmers_in_fastq(
        args.fastq,
        automaton,
        output_count_file=args.output_count,
        output_detail_file=args.output_detail
    )


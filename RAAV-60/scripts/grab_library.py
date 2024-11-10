import gzip

# File paths
barEx_file = '../02_fragBarExt/R1.fastq.gz'
subclean_file = '../01_qualityTrim/subclean.fastq.gz'
output_file = 'matched_sequences.txt'

# Function to read sequences from barEx_file and store in a set
def read_barEx_sequences(file_path):
    sequences = set()
    with gzip.open(file_path, 'rt') as f:
        while True:
            identifier = f.readline().strip()
            if not identifier:
                break
            sequence = f.readline().strip()
            f.readline()  # Skip the '+' line
            f.readline()  # Skip the quality line
            if len(sequence) == 24:
                sequences.add(sequence)
    return sequences

# Read sequences from barEx_file
barEx_sequences = read_barEx_sequences(barEx_file)
print(f"Total barEx sequences read: {len(barEx_sequences)}")

# Open the output file for writing
with open(output_file, 'w') as out_f:
    # Read the subclean file
    with gzip.open(subclean_file, 'rt') as f:
        match_count = 0
        while True:
            identifier = f.readline().strip()
            if not identifier:
                break
            sequence = f.readline().strip()
            f.readline()  # Skip the '+' line
            f.readline()  # Skip the quality line
            # Check if any 24-character sequence from barEx_sequences is in the sequence
            for barEx_seq in barEx_sequences:
                pos = sequence.find(barEx_seq)
                if pos != -1:
                    next_three_chars = sequence[pos + len(barEx_seq):pos + len(barEx_seq) + 3]
                    # Write the matched identifier and the next three characters to the output file
                    out_f.write(f"{identifier}\n{next_three_chars}\n")
                    match_count += 1
                    break
        print(f"Total matches found: {match_count}")

print("Matching sequences have been written to", output_file)
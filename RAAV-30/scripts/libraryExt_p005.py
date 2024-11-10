import gzip

# Ensure the processed_name1 variable is correctly set
processed_name1 = 'AAV23-01_S1_R1_001_LibExt.fastq.gz'
output_file = 'library_005.fastq.gz'

with gzip.open(processed_name1, 'rt') as infile, gzip.open(output_file, 'wt') as outfile:
    while True:
        header = infile.readline().strip()
        if not header:
            break
        sequence = infile.readline().strip()
        plus = infile.readline().strip()
        quality = infile.readline().strip()
        
        # Extract the first three characters from the sequence
        first_three = sequence[:3]
        
        # Write the first three characters to the output file
        outfile.write(f"{header}\n{first_three}\n{plus}\n{quality[:3]}\n")

print(f"First three characters from each sequence saved to {output_file}")
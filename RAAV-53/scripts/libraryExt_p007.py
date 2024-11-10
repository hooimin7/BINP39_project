import sys
import gzip

def extract_first_three_characters(input_file, output_file):
    with gzip.open(input_file, 'rt') as infile, gzip.open(output_file, 'wt') as outfile:
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

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python libraryExt_p007.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    extract_first_three_characters(input_file, output_file)


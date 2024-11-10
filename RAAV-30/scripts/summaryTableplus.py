import pandas as pd
import glob
import os
import gzip
import re


# Define the base directory
base_dir = "/home/hooimin/lu2024-17-19/RAAV-30"

# List of directories ending with _AAV
directories = glob.glob(f"{base_dir}/*_AAV_[0-9]*")

# Initialize an empty DataFrame for the log table
log_table = pd.DataFrame(columns=["Name", "Reads", "Bc", "allBCs", "scBC"])

# Function to count reads in a gzipped FASTQ file
def count_fastq_reads(file_path):
    count = 0
    with gzip.open(file_path, 'rt') as file:
        for line in file:
            count += 1
    # Since each read is represented by 4 lines
    return count // 4

print("log table")
print(log_table)

# Loop through each directory
for dir in directories:
    dir_name = os.path.basename(dir)
    print(f"Processing directory: {dir_name}")

    # Add the directory name to the first column of the log table
    new_row = pd.DataFrame([{"Name": dir_name, "Reads": "NA", "Bc": "NA", "allBCs": "NA", "scBC": "NA"}])
    log_table = pd.concat([log_table, new_row], ignore_index=True)
    index = log_table.index[log_table['Name'] == dir_name].tolist()[0]

    # Define the file pattern for barcodes
    extracted_files = glob.glob(os.path.join(dir, "*_R1_001.fastq.gz"))
    print(f"Found extracted files: {extracted_files}")

    # Count reads using the shell command
    if extracted_files:
        read_count = 0
        for file in extracted_files:
            read_count += count_fastq_reads(file)
        
        log_table.at[index, "Reads"] = read_count
        print(f"Read count for {dir_name}: {read_count}")
    else:
        print(f"No extracted files found for {dir_name}")

    # Define the result file pattern
    result_files = glob.glob(os.path.join(dir, "02_BarExt", "*.out"))
    print(f"Found result files: {result_files}")

    # Extract BCs
    if result_files:
        for file in result_files:
            with open(file) as f:
                for line in f:
                    if "Result:" in line:
                        bc = line.split('\t')[1].strip()
                        log_table.at[index, "Bc"] = bc
                        print(f"BC for {dir_name}: {bc}")
    else:
        print(f"No result files found for {dir_name}")

    # Define the allBCs file pattern
    allBCs_files = glob.glob(os.path.join(dir, "06_traceBack", "02_analysis", "RNA.out"))
    print(f"Found result files: {allBCs_files}")

    # Extract allBCs
    if allBCs_files:
        for file in allBCs_files:
            with open(file) as f:
                for line in f:
                    if "allBCs," in line:
                        match = re.search(r"Original unique barcodes: (\d+)", line)
                        if match:
                            allBCs = match.group(1)
                            log_table.at[index, "allBCs"] = allBCs
                            print(f"allBCs for {dir_name}: {allBCs}")

    # Define the scBC file pattern
    scBC_files = glob.glob(os.path.join(dir, "06_traceBack", "02_analysis", "RNA.out"))
    print(f"Found scBC files: {scBC_files}")
                        
     # Extract scBC
    if scBC_files:
        for file in scBC_files:
            with open(file) as f:
                for line in f:
                    if "scBC," in line:
                        match = re.search(r"SC reduced unique barcodes: (\d+)", line)
                        if match:
                            scBC = match.group(1)
                            log_table.at[index, "scBC"] = scBC
                            print(f"scBC for {dir_name}: {scBC}")                     

# # Save the log table to a CSV file
# log_table.to_csv("01_table/log_table.csv", index=False)

# Save the log table to a CSV file in the specified directory
output_csv_path = "/home/hooimin/lu2024-17-19/RAAV-30/p007_AAV_03/06_traceBack/01_table/log_table.csv"
print(f"Output CSV path: {output_csv_path}")

# Ensure the directory exists
os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
print(f"Directory created: {os.path.dirname(output_csv_path)}")

# Save the DataFrame to the CSV file
log_table.to_csv(output_csv_path, index=False)
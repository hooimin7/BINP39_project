import pandas as pd
import glob
import os


# Define the base directory
base_dir = "/home/hooimin/lu2024-17-19/RAAV-30"

# List of directories ending with _AAV
directories = glob.glob(f"{base_dir}/*_AAV_[0-9]*")

# Initialize an empty DataFrame for the log table
log_table = pd.DataFrame(columns=["Name"])

# Loop through each directory
for dir in directories:
    dir_name = os.path.basename(dir)
    print(f"Processing directory: {dir_name}")

    # Add the directory name to the first column of the log table
    new_row = pd.DataFrame([{"Name": dir_name}])
    log_table = pd.concat([log_table, new_row], ignore_index=True)
    index = log_table.index[log_table['Name'] == dir_name].tolist()[0]


# Save the log table to a CSV file in the specified directory
output_txt_path = "/home/hooimin/lu2024-17-19/RAAV-30/p005_AAV_01/06_traceBack/01_table/loglist_30.txt"
print(f"Output TXT path: {output_txt_path}")

# Ensure the directory exists
os.makedirs(os.path.dirname(output_txt_path), exist_ok=True)
print(f"Directory created: {os.path.dirname(output_txt_path)}")

# Save the DataFrame to the CSV file
log_table.to_csv(output_txt_path, index=False, sep='\t')
import os
import pandas as pd

# Define the directory where the files are located
directory = r"C:\Users\jiang\OneDrive\Research\tc\articles\diazirine\crosslink\LBM\reports-4\LF-type2XL"

# Create a dictionary to store the combined data for each ID
combined_data = {}

# Iterate over all files in the directory
for filename in os.listdir(directory):
    if filename.endswith(".xlsx"):
        filepath = os.path.join(directory, filename)
        file_id = filename.split("_")[-1].split(".")[
            0]  # Extract the ID from the filename

        # Read the file and store its data in the combined_data dictionary
        df = pd.read_excel(filepath)
        if file_id in combined_data:
            combined_data[file_id] = pd.concat([combined_data[file_id], df],
                                               ignore_index=True)
        else:
            combined_data[file_id] = df

output_directory = os.path.join(directory, "combined_data")
os.makedirs(output_directory, exist_ok=True)
# Output the combined data for each ID into separate files
for file_id, df in combined_data.items():
    output_filename = os.path.join(directory + "\combined_data",
                                   file_id + ".xlsx")
    df.to_excel(output_filename, index=False)
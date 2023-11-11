import os
import pandas as pd

# Define Pos Column
column1 = 26  # Column AA, start from 0
column2 = 32  # Column AG, start from 0

# Set the path to the directory containing the Excel files and the directory to save the output files
xlsx_dir = r'C:\Users\jiang\OneDrive\Research\tc\articles\diazirine\crosslink\LBM\reports-4\LF-type2XL\combined_data\duplicates'
output_dir = os.path.join(xlsx_dir, 'duplicates')

# Create the output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Iterate over all files in the directory
for file in os.listdir(xlsx_dir):
    if file.endswith('.xlsx'):
        # Read the Excel file into a DataFrame
        file_path = os.path.join(xlsx_dir, file)
        df = pd.read_excel(file_path)

        # Calculate the repeat times for each row and append it as a new column
        df['Duplicates'] = df.groupby(
            [df.columns[column1],
             df.columns[column2]])[df.columns[0]].transform('count')

        # Drop rows with the same value in the 27th and 33rd columns
        df = df.drop_duplicates(
            subset=[df.columns[column1], df.columns[column2]])

        # Save the modified DataFrame to an Excel file in the output directory
        output_file = os.path.join(output_dir, file)
        df.to_excel(output_file, index=False)
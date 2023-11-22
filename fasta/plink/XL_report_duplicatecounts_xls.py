import os
import pandas as pd

# Define Pos Column
column1 = 24  # start from 0
column2 = 28  # start from 0
rtscolumn = 29  # start from 0

# Set the path to the directory containing the Excel files and the directory to save the output files
xlsx_dir = r'G:\MSdata\202310AHL\Hb\reports\cross-link/'
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
        # Add a new column "allPSMs" and calculate the sum of values in column AM for rows with the same value in column AA and column AG
        df['allPSMs'] = df.groupby([df.columns[column1], df.columns[column2]
                                    ])[df.columns[rtscolumn]].transform('sum')

        # Drop rows with the same value in the 27th and 33rd columns
        df = df.drop_duplicates(
            subset=[df.columns[column1], df.columns[column2]])

        # Save the modified DataFrame to an Excel file in the output directory
        output_file = os.path.join(output_dir, file)
        df.to_excel(output_file, index=False)
import os
import pandas as pd

# Define Pos Column
column1 = 26  # Column AA,start from 0
colume2 = 32  # Column AG,start from 0

# Set the path to the directory containing the CSV files
csv_dir = r'C:\Users\jiang\OneDrive\Research\tc\articles\diazirine\crosslink\BSA\reports-3\BSA-type1XL'

# Iterate over all files in the directory
for file in os.listdir(csv_dir):
    if file.endswith('.csv'):
        # Read the CSV file into a DataFrame
        file_path = os.path.join(csv_dir, file)
        df = pd.read_csv(file_path)

        # Calculate the repeat times for each row and append it as a new column
        df['Repeat Times'] = df.groupby(
            [df.columns[column1],
             df.columns[colume2]])[df.columns[0]].transform('count')
        # Drop rows with the same value in the 27th and 33rd columns
        df = df.drop_duplicates(
            subset=[df.columns[column1], df.columns[colume2]])

        # Save the modified DataFrame to an Excel file
        output_file = os.path.join(csv_dir,
                                   os.path.splitext(file)[0] + '.xlsx')
        df.to_excel(output_file, index=False)
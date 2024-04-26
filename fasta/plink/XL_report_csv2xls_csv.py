import os
import pandas as pd

# Define Pos Column
column1 = 28  # Column AC,start from 0
column2 = 24  # Column Y,start from 0

# Set the path to the directory containing the CSV files
csv_dir = r'G:\MSdata\221027ETF\Trans\reports\cross-link/'

# Iterate over all files in the directory
for file in os.listdir(csv_dir):
    if file.endswith('.csv'):
        # Read the CSV file into a DataFrame
        file_path = os.path.join(csv_dir, file)
        df = pd.read_csv(file_path)

        # Calculate the repeat times for each row and append it as a new column
        df['Repeat Times'] = df.groupby(
            [df.columns[column1],
             df.columns[column2]])[df.columns[0]].transform('count')
        # Drop rows with the same value in the 27th and 33rd columns
        df = df.drop_duplicates(
            subset=[df.columns[column1], df.columns[column2]]
        )  # Cannot delete if subset has reverse value( eg. [111,222] and [222,111]), which not happened in this case (plink alpha pep will always be longer pep)

        # Save the modified DataFrame to an Excel file
        output_file = os.path.join(csv_dir,
                                   os.path.splitext(file)[0] + '.xlsx')
        df.to_excel(output_file, index=False)
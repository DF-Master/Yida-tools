from openpyxl import load_workbook

# Define the file path
file_dir = r'C:\Users\jiang\OneDrive\Research\tc\articles\diazirine\crosslink\LBM\reports-4\LF-type2XL\combined_data\duplicates\report'
file_name = r'residue_counts_dts_1_rts_3.xlsx'
file_path = file_dir + '\\' + file_name
sheet_name = '200.xlsx'

# Load the workbook
workbook = load_workbook(file_path)

# Select the worksheet
worksheet = workbook[sheet_name]

# Define Pos Column !!!!!
column1 = 26  # Column AA,start from 0
column2 = 32  # Column AG,start from 0

# Loop through each row
for row in worksheet.iter_rows(min_row=1, values_only=True):
    # Get the values from Column X and Column AD
    value1 = row[column1]
    value2 = row[column2]

    # Print the formatted output
    print(f"dist /Lactoferrin//A/{value1}/CA, /Lactoferrin//A/{value2}/CA;")
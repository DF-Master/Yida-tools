from Bio import PDB


def cal_distance_pos(animo_pos1,
                     animo_pos2,
                     pdb,
                     name='default',
                     animo_core1='CA',
                     animo_core2='CA',
                     in_one_domain=True,
                     autoprint=False):
    parser = PDB.PDBParser()
    structure = parser.get_structure(name, pdb)
    model = structure[0]
    chain = model['A']
    residue1 = chain[animo_pos1]  # Start from 1 not 0
    residue2 = chain[animo_pos2]
    atoM1 = residue1[animo_core1]
    atom2 = residue2[animo_core2]
    distance = atoM1 - atom2
    if autoprint == True:
        print('Distance between', residue1.get_resname(), animo_core1, 'and',
              residue2.get_resname(), animo_core2, 'is ', distance)
    return distance


import openpyxl
# Define the path to the Excel file
# file_path = r"G:/MSdata/220329ADK/adk/spotlink/data_analyser.xlsx"
file_path = r"G:/MSdata/220329ADK/adk/plink/pLink_task_2022.04.02.20.19.24/reports/data_analyser_plink.xlsx"
pdb_file = r"G:/MSdata/pdb/4ake.pdb"

# Load the workbook
workbook = openpyxl.load_workbook(file_path)

# Select the sheets to compare
sheet_names = workbook.sheetnames

# Colume define
colume1 = 25  # Column Y,start from 1
colume2 = 29  # Column AC

# Distance Column
for sheet_name in sheet_names:
    sheet = workbook[sheet_name]
    max_row = sheet.max_row
    max_col = sheet.max_column

    # Add header for new column
    sheet.cell(row=1, column=max_col + 1, value=pdb_file.split("/")[-1])

    # Iterate through rows to calculate the distance and update the new column
    for row in range(2, max_row + 1):
        pos1 = sheet.cell(row=row, column=colume1).value
        pos2 = sheet.cell(row=row, column=colume2).value
        distance = cal_distance_pos(pos1, pos2, pdb_file)
        sheet.cell(row=row, column=max_col + 1, value=distance)

# Save the modified workbook
workbook.save(file_path)
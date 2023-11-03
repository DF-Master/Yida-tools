import csv

# Specify the input and output file paths
# input_file = r"G:/MSdata/220329ADK/adk/spotlink/WXZ_20220329_AdK_HCDFT_site_filtered.csv"
# output_file = r"G:/MSdata/220329ADK/adk/spotlink/WXZ_20220329_AdK_HCDFT_site_filtered_modified.csv"
input_file = r"G:/MSdata/220329ADK/adk/plink/pLink_task_2022.04.02.20.19.24/reports/adk_2022.04.02.filtered_cross-linked_spectra.csv"
output_file = r"G:/MSdata/220329ADK/adk/plink/pLink_task_2022.04.02.20.19.24/reports/adk_2022.04.02.filtered_cross-linked_spectra_modified.csv"
fasta_file = r"G:/MSdata/fasta/adk.fasta"
split_colume = 4  # start from 0, plink=4, spotlink=2


# Def extract_code(text): (from str(code))
def extract_code(text):
    start_idx = text.index('(') + 1
    end_idx = text.index(')')
    return text[start_idx:end_idx]


# Def read_fasta(fasta_file): (from fasta_file)
def read_fasta(fasta_file):
    with open(fasta_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first row

        # Remove any leading/trailing whitespaces and combine the lines into a single string
        sequence = "".join(line.strip() for line in lines)

    return sequence


# Read the input file and modify the data
with open(input_file, 'r') as file:
    reader = csv.reader(file)
    data = list(reader)  # Convert the reader object into a list
    fasta = read_fasta(fasta_file)

    for row in data:

        text = row[
            split_colume]  # Assuming the third column has the desired string

        if "-" in row[split_colume]:
            # Split the string by "-" into two parts
            str1, str2 = text.split('-')

            # Extract code1 and code2 from the split strings
            code1 = extract_code(str1)
            code2 = extract_code(str2)
            code1_in_fasta = int(code1) + len(
                fasta.split(str1.split('(')[0])[0])
            code2_in_fasta = int(code2) + len(
                fasta.split(str2.split('(')[0])[0])
            # Append the additional columns to the row
            row.extend([
                str1,
                code1,
                str1[int(code1) - 1],
                code1_in_fasta,
                str2,
                code2,
                str2[int(code2) - 1],
                code2_in_fasta,
            ])
        else:
            row.extend([
                "pep1", "pep1(code1)", "pep1[code1]", "code1_in_fasta", "pep2",
                "pep2(code2)", "pep2[code2]", "code2_in_fasta"
            ])

# Write the modified data to the output file
with open(output_file, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(data)
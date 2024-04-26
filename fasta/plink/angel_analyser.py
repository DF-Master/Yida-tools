from Bio.PDB import *


def cal_angel(pdb_file,
              struc_name,
              residue_code1,
              residue_code2,
              residue_code3,
              chain_name='A',
              residue_name1='NZ',
              residue_name2='OH',
              residue_name3='CB'):

    # set up your parser
    parser = PDBParser()
    # get your structure
    structure = parser.get_structure(struc_name, pdb_file)

    # define the model, chain and residue you're interested in
    model = structure[0]
    chain = model[chain_name]  # Replace 'A' with your chain id

    # define three atoms
    res_1 = chain[residue_code1][residue_name1]
    res_2 = chain[residue_code2][residue_name2]
    res_3 = chain[residue_code3][residue_name3]

    # getting coordinates
    vector1 = res_1.get_vector()
    vector2 = res_2.get_vector()
    vector3 = res_3.get_vector()

    #calculate the angle
    angle = calc_angle(vector1, vector2, vector3) * 57.296  # 180/pi = 57.296
    print(angle)
    return angle


if __name__ == '__main__':
    # for i in [293]:
    #     cal_angel("G:/MSdata/pdb/Lactoferrin.pdb",
    #               "Lactoferrin",
    #               i,
    #               181,
    #               181,
    #               chain_name='A',
    #               residue_name1='NZ',
    #               residue_name2='OH',
    #               residue_name3='CB')

    # Open excel file
    import pandas as pd
    file_path = r"G:\MSdata\231215_10MIX\intersection_UNION\cross-link\SUM\intersection_union.xlsx"
    df = pd.read_excel(file_path, sheet_name='T_union')

    # pdb file path
    pdb_dir = 'G:/MSdata/pdb/'
    pdb_dic = {
        'BSA': 'bsa',
        'PPASE': "PPASE",
        'LF': 'Lactoferrin',
        'JYD_CA_': 'conalbumin',
        'ADRM1': 'ADRM1',
        'GFP': 'muGFP',
        'NSP5': 'NSP5',
        'ADK': 'adk',
        'RSV_CA': 'rsv_ca'
    }

    # Save desired columns to lists
    col1 = df.iloc[:, 1].tolist()[1:]
    col23 = df.iloc[:, 23].tolist()[1:]  # Python index starts from 0
    col24 = df.iloc[:, 24].tolist()[1:]
    col27 = df.iloc[:, 27].tolist()[1:]
    col28 = df.iloc[:, 28].tolist()[1:]

    print(col1, col23, col24, col27, col28)

    for i in range(len(col1)):
        for key in pdb_dic:
            if key in col1[i]:
                pdb_name = pdb_dic[key]
                if col27[i] == "K":
                    cal_angel(f"G:/MSdata/pdb/{pdb_name}.pdb",
                              pdb_name,
                              col28[i],
                              col24[i],
                              col24[i],
                              chain_name='A',
                              residue_name1='NZ',
                              residue_name2='OG1',
                              residue_name3='CA')

                elif col23[i] == "K":
                    cal_angel(f"G:/MSdata/pdb/{pdb_name}.pdb",
                              pdb_name,
                              col24[i],
                              col28[i],
                              col28[i],
                              chain_name='A',
                              residue_name1='NZ',
                              residue_name2='OG1',
                              residue_name3='CA')

import csv
from Bio import PDB
import numpy as np
import os
import settings as set


def find_loc_in_fasta(pep, fasta):
    try:
        fasta_pep_split = fasta.split(pep)
    except:
        print('Cannot find pep: ', pep)
    else:
        if len(fasta_pep_split) == 2:
            fasta_length = len(fasta)
            start_pos = len(fasta_pep_split[0]) + 1
            end_pos = len(fasta_pep_split[0]) + len(pep)
            return start_pos, end_pos
        elif len(fasta_pep_split) > 2:
            print("Multi pep Find:", pep)
            print("Please Check:", pep)
            return 1, 1
        elif len(fasta_pep_split) == 1:
            return 1, 1
        else:
            print('Error: ', pep)
            return 1, 1


def find_link_pos(ms_peptide, fasta, spliter="-"):
    pep1, pep2 = ms_peptide.split(spliter)
    pep1_fasta, pep1_num = pep1.strip(')').split('(')
    pep2_fasta, pep2_num = pep2.strip(')').split('(')
    cross_pos_1 = pep1_fasta[int(pep1_num) - 1]
    cross_pos_2 = pep2_fasta[int(pep2_num) - 1]
    pep1_start, pep1_end = find_loc_in_fasta(pep1_fasta, fasta)
    pep2_start, pep2_end = find_loc_in_fasta(pep2_fasta,
                                             fasta)  # if not find return 1,1

    return [[
        pep1_fasta, pep1_num, cross_pos_1, pep1_start, pep1_end,
        pep1_start + int(pep1_num) - 1
    ],
            [
                pep2_fasta, pep2_num, cross_pos_2, pep2_start, pep2_end,
                pep2_start + int(pep2_num) - 1
            ]]


# this function will del repeat element
def cal_repeat_list(list_1, list_2, simple_mode=True):
    print(len(list_1), len(list_2))
    repeat_list = []
    unrepeat_list_1 = []
    unrepeat_list_2 = []

    for i in list_1 + list_2:
        if i in list_2 and i in list_1:
            repeat_list.append(i)
        elif i not in list_2 and i in list_1:
            unrepeat_list_1.append(i)
        elif i in list_2 and i not in list_1:
            unrepeat_list_2.append(i)

    if simple_mode == False:
        try:
            return [
                repeat_list, unrepeat_list_1, unrepeat_list_2,
                len(repeat_list) / (len(unrepeat_list_1) + len(repeat_list)),
                len(repeat_list) / (len(unrepeat_list_2) + len(repeat_list))
            ]
        except:
            return [
                repeat_list, unrepeat_list_1, unrepeat_list_2,
                (len(unrepeat_list_1) + len(repeat_list)),
                (len(unrepeat_list_2) + len(repeat_list))
            ]

    else:
        repeat_list_simple = []
        for i in repeat_list:
            if i not in repeat_list_simple:
                repeat_list_simple.append(i)
        print(len(repeat_list_simple))
        return [repeat_list_simple]


# cal multi repeat list
def cal_multi_list(list_list, simple_mode=True):
    repeat_list = list_list[0]
    for i in range(len(list_list) - 1):
        repeat_list = cal_repeat_list(repeat_list,
                                      list_list[i + 1],
                                      simple_mode=simple_mode)[0]
    return repeat_list


### Analyze amino ratio
# Input list like [["K","L"],["Y","K"],……]
# if pos_form on, input should be [[139,145],[224,269],……]


def report_animo_ratio(
        link_list,
        hard_animo="K",
        pos_count_form=False,
        pos_form=False,
        fasta_file='MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA',
        tableform=False,
        ratio_output=True,
        threshold=0):
    # Trans pos_list to link_list
    if pos_form == True:
        if pos_count_form == True:
            link_list = list(map(list, zip(*link_list)))
            pos_list = link_list[0]
            count_list = link_list[1]
        else:
            pos_list = link_list
            count_list = [1] * len(link_list)
        link_list = []
        for i in pos_list:
            link_list.append(
                [fasta_file[int(i[0]) - 1], fasta_file[int(i[1]) - 1]])
    # print(link_list)
    # Count amino num
    animo_list = list("ACDEFGHIKLMNPQRSTVWY")
    cross_dic = {}
    for animo_tag in animo_list:
        cross_dic[animo_tag] = 0
    cross_dic["all"] = 0
    cross_dic["nofind"] = 0
    n = 0
    for cross_dimer in link_list:
        if threshold == 0:
            threshold = 1008610010
        elif 0 < threshold < 1:
            threshold = threshold * len(link_list)
        n += 1
        if n <= threshold:
            if cross_dimer[0] == hard_animo:
                cross_dic[cross_dimer[1]] += 1 * count_list[n - 1]
                cross_dic["all"] += 1
            elif cross_dimer[1] == hard_animo:
                cross_dic[cross_dimer[0]] += 1 * count_list[n - 1]
                cross_dic["all"] += 1
            else:
                cross_dic["nofind"] += 1
    # Report
    if tableform == False:
        return cross_dic
    # Out dic like {'A': 31, 'C': 26, 'D': 72, 'E': 117, 'F': 11, 'G': 25,……
    else:
        # print('|A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y|ALL|')
        # print('|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|---|')
        tableformdata = []
        if int(cross_dic["all"]) == 0:
            print("No Crosslink")
        else:
            for i in cross_dic:
                if i == 'nofind':
                    continue
                elif i != 'all' and i != 'nofind':
                    if ratio_output == True:
                        tableformdata.append(
                            format(cross_dic[i] / cross_dic['all'], '.3f'))
                    else:
                        tableformdata.append(format(cross_dic[i], '.1f'))
                elif i == 'all':
                    tableformdata.append(str(cross_dic['all']) + "|")

            tableformdata = '|'.join(tableformdata)
            print('|' + tableformdata)
        return cross_dic


### Plink analyze
# Unit Form


def plink2normalform(plink_file, output_dir, del_mode=True):
    folder = os.path.exists(output_dir)

    if not folder:
        os.makedirs(output_dir)  # makedirs 创建文件时如果路径不存在会创建这个路径
        print('文件夹创建成功：', output_dir)

    else:
        print('文件夹已经存在：', output_dir)

    with open(plink_file) as f:
        reader = csv.reader(f)
        for line in reader:
            open(output_dir + str(line[1].split(".")[0]) + ".csv",
                 "a",
                 newline='')
            with open(output_dir + str(line[1].split(".")[0]) + ".csv",
                      "r") as outfile:
                csv_reader_f = csv.reader(outfile)
                if not any(csv_reader_f):
                    with open(output_dir + str(line[1].split(".")[0]) + ".csv",
                              "a",
                              newline='') as outfile_append:
                        writer = csv.writer(outfile_append)
                        if del_mode == False:
                            writer.writerow([
                                "Order", "Title", "Charge", "Precursor_Mass",
                                "Peptide", "Peptide_Type", "Linker",
                                "Peptide_Mass", "Modifications", "Evalue",
                                "Score", "Precursor_Mass_Error(Da)",
                                "Precursor_Mass_Error(ppm)", "Proteins",
                                "Protein_Type", "FileID", "LabelID",
                                "Alpha_Matched", "Beta_Matched",
                                "Alpha_Evalue", "Beta_Evalue"
                            ])
                else:
                    with open(output_dir + str(line[1].split(".")[0]) + ".csv",
                              "a",
                              newline='') as outfile_append:
                        writer = csv.writer(outfile_append)
                        if "." in line[1]:
                            if del_mode == True:
                                writer.writerow([
                                    line[0], line[1], line[2], line[3],
                                    line[5], line[4], line[7], line[8],
                                    line[13], line[14], line[9], line[10]
                                ])
                            else:
                                writer.writerow(line)


#  Input ms_csv_file
def plink_report_link_pos(ms_csv_file,
                          fasta,
                          ms_peptide_list=5,
                          type="crosslink"):
    link_pos_list = []
    with open(ms_csv_file, "r") as f:
        reader = csv.reader(f)
        for target in reader:
            if len(target) >= 2:

                try:
                    if type == "crosslink":
                        link_pos_list.append(
                            find_link_pos(target[ms_peptide_list], fasta))
                    elif type == "monolink":
                        link_pos_list.append(
                            find_link_pos(
                                target[ms_peptide_list] + "-" +
                                str(fasta[2:-2]) + "(1)", fasta))
                    elif type == "looplink":
                        link_pos_list.append(
                            find_link_pos(
                                str(")-" +
                                    target[ms_peptide_list].split("(")[0] +
                                    "(").join(
                                        target[ms_peptide_list].split(")(")),
                                fasta))
                    else:
                        print("Wrong Type")
                        break

                    link_pos_list[-1].append([
                        target[ms_peptide_list + 5],
                        target[ms_peptide_list + 6]
                    ])
                    link_pos_list = [
                        i for i in sorted(link_pos_list,
                                          key=lambda x: x[-1][-1],
                                          reverse=True)
                        # if float(i[-1][-1]) <= 0.25
                    ]
                except:
                    if target[
                            ms_peptide_list] != "Alpha Peptide Protein" and target[
                                ms_peptide_list] != "Peptide":
                        print("Target link position cannot be analyzed:",
                              target[ms_peptide_list])
            else:
                continue
    return link_pos_list


def plink_report_valid_link_sfDR(
        ms_csv_dir,
        fasta="MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA",
        filtermode=True,
        least_ssn=2,
        scorek=0.05,
        threshold=1,
        ms_peptide_list=5,
        type="crosslink",
        score_count_mode=False):

    report_list = plink_report_link_pos(ms_csv_dir,
                                        fasta,
                                        ms_peptide_list=ms_peptide_list,
                                        type=type)

    if filtermode == True:
        report_list = [
            i
            for i in sorted(report_list, key=lambda x: x[2][1], reverse=False)
            [0:int(threshold * len(report_list))]
        ]

    report_list = [[max([i[0][-1], i[1][-1]]),
                    min([i[0][-1], i[1][-1]])] for i in report_list]

    count_loc_list = [[], []]
    for i in report_list:
        if i not in count_loc_list[0]:
            count_loc_list[0].append(i)
            count_loc_list[1].append(1)
        else:
            count_loc_list[1][count_loc_list[0].index(i)] += 1
    count_loc_list = list(map(list, zip(*count_loc_list)))
    if score_count_mode == True:
        for i in count_loc_list:
            i[1] = 1 - 3**(-scorek * i[1])
    hard_loc_list = []
    if least_ssn > 0:
        for i in count_loc_list:
            if i[1] >= least_ssn:
                hard_loc_list.append(i[0])
    return [hard_loc_list, report_list, count_loc_list]


#####
# Distance
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


def cal_distance_pos_list(
        pos_list,
        fasta="MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA",
        pdb_file='G:/MSdata/pdb/bsa.pdb',
        animo_core1='CA',
        animo_core2='CA',
        max_distance=20,
        min_distance=6,
        threshold=1,
        print_all=True,
        read_mod=False):

    if threshold > 1:
        pos_list = pos_list[0:threshold]
    elif 0 < threshold < 1:
        pos_list = pos_list[0:int(threshold * len(pos_list))]

    dis_list = [[], []]
    over_max_num = 0
    over_min_num = 0
    for i in pos_list:
        dis = float(
            format(
                cal_distance_pos(int(i[0]),
                                 int(i[1]),
                                 pdb_file,
                                 animo_core1=animo_core1,
                                 animo_core2=animo_core2), '.5f'))
        if dis >= max_distance:
            over_max_num += 1
        elif dis <= min_distance:
            over_min_num += 1
        if print_all == True:
            print(i, dis, fasta[int(i[0]) - 1], fasta[int(i[1]) - 1])
        dis_list[0].append(i)
        dis_list[1].append(dis)

    if read_mod == True:
        print("| Average |", format(np.mean(dis_list[1]), '.5f'), '| Median |',
              format(np.median(dis_list[1]),
                     '.5f'), '| Std |', np.std(dis_list[1]), '| Num |',
              len(dis_list[1]), '| Over Max Num |', over_max_num,
              '| Over Min Num |', over_min_num, "|")
        dis_list = list(map(list, zip(*dis_list)))

        return dis_list
    else:
        if print_all == True:
            print(dis_list[1])
        return [
            dis_list[0], dis_list[1], "| Average |",
            format(np.mean(dis_list[1]), '.5f'), '| Median |',
            format(np.median(dis_list[1]), '.5f'), '| Std |',
            np.std(dis_list[1]), '| Num |',
            len(dis_list[1]), '| Over Max Num |', over_max_num,
            '| Over Min Num |', over_min_num, "|"
        ]


#####
# Fasta analyze
def calc_fasta(fasta, table_form=False):
    animo_list = list("ACDEFGHIKLMNPQRSTVWY")
    cross_dic = {}

    for animo_tag in animo_list:
        cross_dic[animo_tag] = 0
    cross_dic["all"] = 0
    cross_dic["nofind"] = 0

    for i in list(fasta):
        if i in cross_dic:
            cross_dic[i] += 1
            cross_dic['all'] += 1
        else:
            cross_dic['nofind'] += 1
            cross_dic['all'] += 1

    if table_form == True:
        # print('|A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y|ALL|')
        # print('|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|---|')
        tableformdata = []
        for i in cross_dic:
            if i != 'all':
                tableformdata.append(
                    format(cross_dic[i] / cross_dic['all'], '.3f'))
            elif i == 'all':
                tableformdata.append(str(cross_dic['all']) + "|")
                break
        tableformdata = '|'.join(tableformdata)
        print('|' + tableformdata)
        return cross_dic
    else:
        return cross_dic


# Neighborhood analyze
def analyze_neighborhood(fasta,
                         target='K',
                         distance=7,
                         percent=False,
                         table_form=False):
    animo_list = list("ACDEFGHIKLMNPQRSTVWY")
    cross_dic = {}
    for animo_tag in animo_list:
        cross_dic[animo_tag] = 0
    cross_dic["all"] = 0
    cross_dic["nofind"] = 0

    fasta_list = list(fasta)
    for i in range(len(fasta)):
        if fasta[i] == target:
            for j in range(i - distance, i + distance + 1):
                try:
                    if fasta_list[j] in animo_list:
                        cross_dic[fasta_list[j]] += 1
                        cross_dic["all"] += 1
                    else:
                        cross_dic['nofind'] += 1

                except:
                    continue
    if table_form == True:
        # print('|A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y|ALL|')
        # print('|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|---|')
        tableformdata = []
        for i in cross_dic:
            if i != 'all':
                tableformdata.append(
                    format(cross_dic[i] / cross_dic['all'], '.3f'))
            elif i == 'all':
                tableformdata.append(str(cross_dic['all']) + "|")
                break
        tableformdata = '|'.join(tableformdata)
        print('|' + tableformdata)
        return cross_dic

    elif percent == False:
        return (cross_dic)
    else:
        for animo_tag in animo_list:
            cross_dic[animo_tag] = cross_dic[animo_tag] / cross_dic['all']
        return (cross_dic)


def analyze_distance_neighboorhood(fasta,
                                   pdb_file,
                                   target='K',
                                   distance_max=25,
                                   distance_min=6):
    animo_list = list("ACDEFGHIKLMNPQRSTVWY")
    dis_neighboorhood_list = [[i, []] for i in animo_list]
    fasta_list = list(fasta)
    num = 0
    for i in fasta_list:
        num += 1
        if i == target:
            for j in range(len(fasta_list)):
                dis_i_j = float(
                    format(cal_distance_pos(num, j + 1, pdb_file), '.5f'))
                if dis_i_j <= distance_max and dis_i_j >= distance_min:
                    dis_neighboorhood_list[animo_list.index(
                        fasta_list[j])][1].append([[num, j + 1], dis_i_j])
    for i in dis_neighboorhood_list:
        try:
            dis_list = list(map(list, zip(*i[1])))
            i.append([
                str(i[0]), "Average",
                format(np.mean(dis_list[1]), '.5f'), ', Median',
                format(np.median(dis_list[1]), '.5f'), ', Std',
                np.std(dis_list[1]), ', Num',
                len(dis_list[1])
            ])
            print(i[2])
        except:
            print("Error")

    return (dis_neighboorhood_list)


if __name__ == "__main__":
    plink2normalform(
        "G:/MSdata/230703BQMCOLD/LBM/reports/Lactoferrin_con_2023.07.03.filtered_cross-linked_spectra.csv",
        "G:/MSdata/230703BQMCOLD/LBM/reports/crosslink_withdis/",
        del_mode=False)
    # plink2normalform(
    #     "G:/MSdata/230703BQMCOLD/BSA/reports/bsa_con_2023.07.03.filtered_cross-linked_spectra.csv",
    #     "G:/MSdata/230703BQMCOLD/BSA/reports/crosslink_withdis/",
    #     del_mode=False)
    # plink2normalform(
    #     "G:/MSdata/230703BQMCOLD/CA/reports/conalbumin_con_2023.07.04.filtered_cross-linked_spectra.csv",
    #     "G:/MSdata/230703BQMCOLD/CA/reports/crosslink_withdis/",
    #     del_mode=False)
    # plink2normalform(
    #     "C:/Users/jiang/OneDrive/Research/tc/articles/diazirine/crosslink/BQ/reports-4/conalbumin_con_2023.03.24.filtered_loop-linked_spectra.csv",
    #     "C:/Users/jiang/OneDrive/Research/tc/articles/diazirine/crosslink/BQ/reports-4/looplink_withdis/",
    #     del_mode=False)
    # plink2normalform(
    #     "C:/Users/jiang/OneDrive/Research/tc/articles/diazirine/crosslink/LBM/reports-4/Lactoferrin_con_2023.03.23.filtered_cross-linked_spectra.csv",
    #     "C:/Users/jiang/OneDrive/Research/tc/articles/diazirine/crosslink/LBM/reports-4/crosslink_withdis/",
    #     del_mode=False)
    # plink2normalform(
    #     "C:/Users/jiang/OneDrive/Research/tc/articles/diazirine/crosslink/LBM/reports-4/Lactoferrin_con_2023.03.23.filtered_loop-linked_spectra.csv",
    #     "C:/Users/jiang/OneDrive/Research/tc/articles/diazirine/crosslink/LBM/reports-4/looplink_withdis/",
    #     del_mode=False)

    # with open(set.fasta_root + "Lactoferrin" + ".fasta", "r") as f:
    #     fasta = "".join([i.strip() for i in f][1:])
    # print(fasta)
    # analyze_distance_neighboorhood(
    #     fasta, set.pdb_root + "LF-AF-Q6LBN7-F1-model_v4.pdb")

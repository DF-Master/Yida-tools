import csv
from sys import float_repr_style
from Bio import PDB
import numpy as np

#######
# Read ms_data and analyse crosslink information


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
        else:
            # print("Find pep Fail:", pep)
            return 1, 1


def find_link_pos(ms_peptide, fasta):
    pep1, pep2 = ms_peptide.split('-')
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

    for i in list_1:
        if i in list_2 and i not in repeat_list:
            repeat_list.append(i)
        elif i not in list_2 and i not in repeat_list:
            unrepeat_list_1.append(i)
    unrepeat_list_2 = [
        i for i in list_2 if i not in unrepeat_list_2 and i not in repeat_list
    ]

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
        try:
            print(
                len(repeat_list) / (len(unrepeat_list_1) + len(repeat_list)),
                len(repeat_list) / (len(unrepeat_list_2) + len(repeat_list)),
                len(repeat_list), (len(unrepeat_list_1) + len(repeat_list)),
                (len(unrepeat_list_2) + len(repeat_list)))
            return [repeat_list]
        except:
            print((len(unrepeat_list_1) + len(repeat_list)),
                  (len(unrepeat_list_2) + len(repeat_list)), len(repeat_list),
                  (len(unrepeat_list_1) + len(repeat_list)),
                  (len(unrepeat_list_2) + len(repeat_list)))
            return [repeat_list]


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
        pos_form=False,
        fasta_file='MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA',
        tableform=False,
        ratio_output=True,
        threshold=0):
    # Trans pos_list to link_list
    if pos_form == True:
        pos_list = link_list
        link_list = []
        for i in pos_list:
            link_list.append(
                [fasta_file[int(i[0]) - 1], fasta_file[int(i[1]) - 1]])
    # Count amino num
    animo_list = list("ACDEFGHIKLMNPQRSTVWY")
    cross_dic = {}
    for animo_tag in animo_list:
        cross_dic[animo_tag] = 0
    cross_dic["all"] = 0
    cross_dic["nofind"] = 0
    if threshold == 0:
        for cross_dimer in link_list:
            if cross_dimer[0] == hard_animo:
                cross_dic[cross_dimer[1]] += 1
                cross_dic["all"] += 1
            elif cross_dimer[1] == hard_animo:
                cross_dic[cross_dimer[0]] += 1
                cross_dic["all"] += 1
            else:
                cross_dic["nofind"] += 1
    elif threshold > 1:
        n = 0
        for cross_dimer in link_list:
            n += 1
            if n <= threshold:
                if cross_dimer[0] == hard_animo:
                    cross_dic[cross_dimer[1]] += 1
                    cross_dic["all"] += 1
                elif cross_dimer[1] == hard_animo:
                    cross_dic[cross_dimer[0]] += 1
                    cross_dic["all"] += 1
                else:
                    cross_dic["nofind"] += 1

            else:
                break
    elif 0 < threshold < 1:
        n = 0
        for cross_dimer in link_list:
            n += 1
            if n <= threshold * len(link_list):
                if cross_dimer[0] == hard_animo:
                    cross_dic[cross_dimer[1]] += 1
                    cross_dic["all"] += 1
                elif cross_dimer[1] == hard_animo:
                    cross_dic[cross_dimer[0]] += 1
                    cross_dic["all"] += 1
                else:
                    cross_dic["nofind"] += 1

            else:
                break

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


def report_score(
        spotlink_report_valid_link_sfDR_output,
        hard_animo="K",
        pos_form=False,
        fasta_file='MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA',
        tableform=False,
        ratio_output=True):

    link_list = spotlink_report_valid_link_sfDR_output[0]
    report_list = spotlink_report_valid_link_sfDR_output[1]

    if pos_form == True:
        pos_list = link_list
        link_list = []
        for i in pos_list:
            link_list.append(
                [fasta_file[int(i[0]) - 1], fasta_file[int(i[1]) - 1]])

    # Count amino num
    animo_list = list("ACDEFGHIKLMNPQRSTVWY")
    cross_dic = {}
    for animo_tag in animo_list:
        cross_dic[animo_tag] = float(0)
    cross_dic["all"] = 0
    cross_dic["nofind"] = 0
    for cross_dimer in link_list:
        if cross_dimer[0] == hard_animo:

            cross_dic[cross_dimer[1]] += float(
                report_list[cross_dic["all"]][2][1])
            cross_dic["all"] += 1
        elif cross_dimer[1] == hard_animo:
            cross_dic[cross_dimer[0]] += float(
                report_list[cross_dic["all"]][2][1])
            cross_dic["all"] += 1
        else:
            cross_dic["nofind"] += 1
    # Report
    print(cross_dic)

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
                        tableformdata.append(format(cross_dic[i], '.3f'))
                    else:
                        tableformdata.append(format(cross_dic[i], '.1f'))
                elif i == 'all':
                    tableformdata.append(str(cross_dic['all']) + "|")

            tableformdata = '|'.join(tableformdata)
            print('|' + tableformdata)
        return cross_dic


#####
# Spotlink analyze


def fdr2sfdr(file_name):
    # Spotlink
    return 'site'.join(file_name.split('result'))


# Input ms_csv_file
def spotlink_report_link_pos(ms_csv_file, fasta, ms_peptide_list=5, f2sf=True):
    # Spotlink
    link_pos_list = []
    with open(ms_csv_file, "r") as f:
        reader = csv.reader(f)
        for target in reader:
            try:
                link_pos_list.append(
                    find_link_pos(target[ms_peptide_list], fasta))
                if f2sf == True:
                    link_pos_list[-1].append([
                        target[ms_peptide_list + 5],
                        target[ms_peptide_list + 6],
                        target[ms_peptide_list + 7]
                    ])
                else:
                    link_pos_list[-1].append([
                        target[ms_peptide_list + 7],
                        target[ms_peptide_list + 8],
                        target[ms_peptide_list + 9],
                    ])
            except:
                if target[
                        ms_peptide_list] != "Alpha Peptide Protein" and target[
                            ms_peptide_list] != "Peptide":
                    print("Target link position cannot be analyzed:",
                          target[ms_peptide_list])
    return link_pos_list


def spotlink_report_valid_link_sfDR(
        ms_csv_dir,
        fasta="MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA",
        set=False,
        filtermode=True,
        least_ssn=1,
        threshold=800,
        ms_peptide_list=5,
        f2sf=False):
    # Spotlink
    if f2sf == True:
        ms_csv_dir = fdr2sfdr(ms_csv_dir)
        ms_peptide_list = 2

    report_list = spotlink_report_link_pos(ms_csv_dir,
                                           fasta,
                                           ms_peptide_list=ms_peptide_list,
                                           f2sf=f2sf)
    if threshold < 1:
        tsd = int(threshold * len(report_list))
    else:
        tsd = threshold
    if filtermode == True and f2sf == True:
        report_list = [
            i for i in sorted(
                report_list, key=lambda x: x[2][1], reverse=False)[0:tsd]
            if int(i[2][2]) >= least_ssn and float(i[2][1]) <= 0.01
        ]
    elif filtermode == True and f2sf == False:
        report_list = [
            i for i in sorted(report_list, key=lambda x: x[2][1], reverse=True)
            [0:tsd] if float(i[2][2]) <= 0.01
        ]
    #print(report_list)
    if set == False:
        hard_loc_list = [[
            max([i[0][-1], i[1][-1]]),
            min([i[0][-1], i[1][-1]])
        ] for i in report_list]
    elif set == True:
        hard_loc_list_all = [[
            max([i[0][-1], i[1][-1]]),
            min([i[0][-1], i[1][-1]])
        ] for i in report_list]
        hard_loc_list = []
        for i in hard_loc_list_all:
            if i not in hard_loc_list:
                hard_loc_list.append(i)

    # report_score([hard_loc_list, report_list], pos_form=True, tableform=True)

    # report_animo_ratio(hard_loc_list,
    #                    pos_form=True,
    #                    fasta_file=fasta,
    #                    tableform=True)

    return [hard_loc_list, report_list]


### Plink analyze
# Unit Form


def plink2crosslink(plink_file, output_dir):
    import os
    folder = os.path.exists(output_dir)
    
    if not folder:
        os.makedirs(output_dir)  # makedirs 创建文件时如果路径不存在会创建这个路径
        print('文件夹创建成功：', output_dir)

    else:
        print('文件夹已经存在：', output_dir)
    
    with open(plink_file) as f:
        reader = csv.reader(f)
        for line in reader:
            if "." in line[1]:
                with open(output_dir + str(line[1].split(".")[0]) + ".csv",
                          "a",
                          newline='') as outfile:
                    writer = csv.writer(outfile)
                    writer.writerow([
                        line[0], line[1], line[2], line[3], line[5], line[4],
                        line[7], line[8], line[13], line[14], line[9], line[10]
                    ])


# Input ms_csv_file
def plink_report_link_pos(ms_csv_file, fasta, ms_peptide_list=5):
    # Spotlink
    link_pos_list = []
    with open(ms_csv_file, "r") as f:
        reader = csv.reader(f)
        for target in reader:
            if len(target) >= 2:

                try:
                    link_pos_list.append(
                        find_link_pos(target[ms_peptide_list], fasta))

                    link_pos_list[-1].append([
                        target[ms_peptide_list + 5],
                        target[ms_peptide_list + 6]
                    ])
                    link_pos_list = [
                        i for i in sorted(link_pos_list,
                                          key=lambda x: x[-1][-1],
                                          reverse=True)
                        if float(i[-1][-1]) <= 0.5
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
        set=True,
        filtermode=True,
        least_ssn=2,
        threshold=1,
        ms_peptide_list=5):

    report_list = plink_report_link_pos(ms_csv_dir,
                                        fasta,
                                        ms_peptide_list=ms_peptide_list)

    if filtermode == True:
        report_list = [
            i
            for i in sorted(report_list, key=lambda x: x[2][1], reverse=False)
            [0:int(threshold * len(report_list))]
        ]
        #print(report_list)
        # report_list = [
        #     i for i in sorted(report_list, key=lambda x: x[2][1], reverse=True)
        #     [0:100]
        # ]

    if set == False:
        hard_loc_list = [[
            max([i[0][-1], i[1][-1]]),
            min([i[0][-1], i[1][-1]])
        ] for i in report_list]
    elif set == True:
        hard_loc_list_all = [[
            max([i[0][-1], i[1][-1]]),
            min([i[0][-1], i[1][-1]])
        ] for i in report_list]

    count_loc_list = [[], []]
    for i in hard_loc_list_all:
        if i not in count_loc_list[0]:
            count_loc_list[0].append(i)
            count_loc_list[1].append(1)
        else:
            count_loc_list[1][count_loc_list[0].index(i)] += 1
    count_loc_list = list(map(list, zip(*count_loc_list)))

    hard_loc_list = []
    for i in count_loc_list:
        if i[1] >= least_ssn:
            hard_loc_list.append(i[0])
    # report_animo_ratio(hard_loc_list,
    #                    pos_form=True,
    #                    fasta_file=fasta,
    #                    tableform=True)
    return [hard_loc_list, report_list]


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
            format(cal_distance_pos(int(i[0]), int(i[1]), pdb_file), '.5f'))
        if dis >= max_distance:
            over_max_num += 1
        elif dis <= min_distance:
            over_min_num += 1
        if print_all == True:
            print(i, dis, fasta[int(i[0]) - 1], fasta[int(i[1]) - 1])
        dis_list[0].append(i)
        dis_list[1].append(dis)

    print("| Average |", format(np.mean(dis_list[1]), '.5f'), '| Median |',
          format(np.median(dis_list[1]),
                 '.5f'), '| Std |', np.std(dis_list[1]), '| Num |',
          len(dis_list[1]), '| Over Max Num |', over_max_num,
          '| Over Min Num |', over_min_num, "|")
    if read_mod == True:
        dis_list = list(map(list, zip(*dis_list)))

        return dis_list
    else:
        if print_all == True:
            print(dis_list[1])
        return dis_list[1]


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
                                   distance_max=20,
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


def cal_distance(animo_pos1,
                 animo_pos2,
                 pdb,
                 name='default',
                 animo_core1='CA',
                 animo_core2='CA',
                 in_one_domain=True,
                 autoprint=True):
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


def mod_fasta_head(fasta_str):
    # print(len(fasta_str))
    fasta_str = "K" + fasta_str[1:-1] + "K"
    # print(len(fasta_str))
    return fasta_str


if __name__ == '__main__':
    plink2crosslink(
        "G:/MSdata/230303BQM/BSA/reports/bsa_con_2023.03.03.filtered_mono-linked_spectra.csv",
        "G:/MSdata/230303BQM/BSA/reports/monolink_spotlinkform/")
    plink2crosslink(
        "G:/MSdata/230303BQM/BQ/reports/conalbumin_con_2023.03.03.filtered_mono-linked_spectra.csv",
        "G:/MSdata/230303BQM/BQ/reports/monolink_spotlinkform/")
    plink2crosslink(
        "G:/MSdata/230303BQM/LBM/reports/Lactoferrin_con_2023.03.03.filtered_mono-linked_spectra.csv",
        "G:/MSdata/230303BQM/LBM/reports/monolink_spotlinkform/")
    
    
    # print("*" * 40)
    # print("--BSA--" * 5)
    # print("*" * 40)
    # ## define part
    # plink_root_dir = "G:/MSdata/230303BQM/BSA/reports/crosslink_spotlinkform/"
    # ### fasta
    # bsa_fasta = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA"
    # ### file_dir

    # plink_file_230112_B0_fDR = plink_root_dir + "JYD_20230112_BSA0.csv"
    # plink_file_230112_B1_fDR = plink_root_dir + "JYD_20230112_BSA1.csv"
    # plink_file_230112_B2_fDR = plink_root_dir + "JYD_20230112_BSA2.csv"
    # plink_file_230112_B5_fDR = plink_root_dir + "JYD_20230112_BSA5.csv"
    # plink_file_230112_B10_fDR = plink_root_dir + "JYD_20230112_BSA10.csv"
    # plink_file_230112_B25_fDR = plink_root_dir + "JYD_20230112_BSA25.csv"
    # plink_file_230112_B50_fDR = plink_root_dir + "JYD_20230112_BSA50.csv"
    # plink_file_230112_B100_fDR = plink_root_dir + "JYD_20230112_BSA100.csv"
    # plink_file_230112_B200_fDR = plink_root_dir + "JYD_20230112_BSA200.csv"

    # plink_root_dir = "G:/MSdata/220913BSABQXL/plink/pLink_task_2022.09.15.14.35.11/reports/crosslink_spotlinkform/"

    # plink_file_0908_B0_fDR = plink_root_dir + "JYD_20220908_B0.csv"
    # plink_file_0908_B5_fDR = plink_root_dir + "JYD_20220908_B5.csv"
    # plink_file_0908_B25_fDR = plink_root_dir + "JYD_20220908_B25.csv"
    # plink_file_0908_B100_fDR = plink_root_dir + "JYD_20220908_B100.csv"
    # plink_file_0908_B1_fDR = plink_root_dir + "JYD_20220908_B1.csv"

    # B0_0908_p = plink_report_valid_link_sfDR(plink_file_0908_B0_fDR,
    #                                          bsa_fasta)[0]
    # B1_0908_p = plink_report_valid_link_sfDR(plink_file_0908_B1_fDR,
    #                                          bsa_fasta)[0]
    # B5_0908_p = plink_report_valid_link_sfDR(plink_file_0908_B5_fDR,
    #                                          bsa_fasta)[0]
    # B25_0908_p = plink_report_valid_link_sfDR(plink_file_0908_B25_fDR,
    #                                           bsa_fasta)[0]
    # B100_0908_p = plink_report_valid_link_sfDR(plink_file_0908_B100_fDR,
    #                                            bsa_fasta)[0]

    # B0_230112_p = plink_report_valid_link_sfDR(plink_file_230112_B0_fDR,
    #                                          bsa_fasta)[0]
    # B1_230112_p = plink_report_valid_link_sfDR(plink_file_230112_B1_fDR,
    #                                          bsa_fasta)[0]
    # B2_230112_p = plink_report_valid_link_sfDR(plink_file_230112_B2_fDR,
    #                                          bsa_fasta)[0]
    # B5_230112_p = plink_report_valid_link_sfDR(plink_file_230112_B5_fDR,
    #                                          bsa_fasta)[0]
    # B10_230112_p = plink_report_valid_link_sfDR(plink_file_230112_B10_fDR,
    #                                          bsa_fasta)[0]
    # B25_230112_p = plink_report_valid_link_sfDR(plink_file_230112_B25_fDR,
    #                                           bsa_fasta)[0]
    # B50_230112_p = plink_report_valid_link_sfDR(plink_file_230112_B50_fDR,
    #                                          bsa_fasta)[0]
    # B100_230112_p = plink_report_valid_link_sfDR(plink_file_230112_B100_fDR,
    #                                            bsa_fasta)[0]
    # B200_230112_p = plink_report_valid_link_sfDR(plink_file_230112_B200_fDR,
    #                                            bsa_fasta)[0]
    # ## repeat part

    # # B0_repeat = cal_multi_list([B0_230112_p, B0_230112_p])
    # # B1_repeat = cal_multi_list([B1_230112_p, B1_230112_p])
    # # B2_repeat = cal_multi_list([B2_230112_p, B2_230112_p])
    # # B5_repeat = cal_multi_list([B5_230112_p, B5_230112_p])
    # # B10_repeat = cal_multi_list([B10_230112_p, B10_230112_p])
    # # B25_repeat = cal_multi_list([B25_230112_p, B25_230112_p])
    # # B50_repeat = cal_multi_list([B50_230112_p, B50_230112_p])
    # # B100_repeat = cal_multi_list([B100_230112_p, B100_230112_p])
    # # B200_repeat = cal_multi_list([B200_230112_p, B200_230112_p])



    # B0_repeat = cal_multi_list([B0_230112_p, B0_0908_p])
    # B1_repeat = cal_multi_list([B1_230112_p, B1_0908_p])
    # B5_repeat = cal_multi_list([B5_230112_p, B5_0908_p])
    # B25_repeat = cal_multi_list([B25_230112_p, B25_0908_p])
    # B100_repeat = cal_multi_list([B100_230112_p, B100_0908_p])

    # # B_repeat_list = [B0_repeat, B1_repeat, B2_repeat,B5_repeat,B10_repeat, B25_repeat, B50_repeat,B100_repeat,B200_repeat]
    # B_repeat_list = [B0_repeat, B1_repeat, B5_repeat, B25_repeat, B100_repeat]

    # print('|-A-|-C-|-D-|-E-|-F-|-G-|-H-|-I-|-K-|-L-|-M-|-N-|-P-|-Q-|-R-|-S-|-T-|-V-|-W-|-Y-|ALL|')
    # print('|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|---|')
    # for i in B_repeat_list:
    #     k = [
    #         j for j in i
    #         if bsa_fasta[int(j[0]) - 1] == "K" or bsa_fasta[int(j[1]) -
    #                                                         1] == "K"
    #     ]
    #     # cal_distance_pos_list(k,
    #     #                       fasta=bsa_fasta,
    #     #                       pdb_file='G:/MSdata/pdb/bsa.pdb',
    #     #                       print_all=False,
    #     #                       threshold=60)
    #     report_animo_ratio(k,
    #                        pos_form=True,
    #                        fasta_file=bsa_fasta,
    #                        tableform=True,
    #                        ratio_output=False,
    #                        threshold=600)
    # for i in B_repeat_list:
    #     k = [
    #         j for j in i
    #         if bsa_fasta[int(j[0]) - 1] == "K" or bsa_fasta[int(j[1]) -
    #                                                         1] == "K"
    #     ]
    #     cal_distance_pos_list(k,
    #                           fasta=bsa_fasta,
    #                           pdb_file='G:/MSdata/pdb/bsa.pdb',
    #                           print_all=False,
    #                           threshold=60)
    #     # report_animo_ratio(k,
    #     #                    pos_form=True,
    #     #                    fasta_file=bsa_fasta,
    #     #                    tableform=True,
    #     #                    ratio_output=False,
    #     #                    threshold=600)


    # print("*" * 40)
    # print("--conalbumin--" * 5)
    # print("*" * 40)
    # ## define part
    # plink_root_dir = "G:/MSdata/230303BQM/BQ/reports/crosslink_spotlinkform/"
    # ### fasta
    # bq_fasta = "MKLILCTVLSLGIAAVCFAAPPKSVIRWCTISSPEEKKCNNLRDLTQQERISLTCVQKATYLDCIKAIANNEADAISLDGGQAFEAGLAPYKLKPIAAEVYEHTEGSTTSYYAVAVVKKGTEFTVNDLQGKTSCHTGLGRSAGWNIPIGTLLHRGAIEWEGIESGSVEQAVAKFFSASCVPGATIEQKLCRQCKGDPKTKCARNAPYSGYSGAFHCLKDGKGDVAFVKHTTVNENAPDQKDEYELLCLDGSRQPVDNYKTCNWARVAAHAVVARDDNKVEDIWSFLSKAQSDFGVDTKSDFHLFGPPGKKDPVLKDLLFKDSAIMLKRVPSLMDSQLYLGFEYYSAIQSMRKDQLTPSPRENRIQWCAVGKDEKSKCDRWSVVSNGDVECTVVDETKDCIIKIMKGEADAVALDGGLVYTAGVCGLVPVMAERYDDESQCSKTDERPASYFAVAVARKDSNVNWNNLKGKKSCHTAVGRTAGWVIPMGLIHNRTGTCNFDEYFSEGCAPGSPPNSRLCQLCQGSGGIPPEKCVASSHEKYFGYTGALRCLVEKGDVAFIQHSTVEENTGGKNKADWAKNLQMDDFELLCTDGRRANVMDYRECNLAEVPTHAVVVRPEKANKIRDLLERQEKRFGVNGSEKSKFMMFESQNKDLLFKDLTKCLFKVREGTTYKEFLGDKFYTVISSLKTCNPSDILQMCSFLEGK"
    # ### file_dir
    # plink_file_230206_BQ1_fDR = plink_root_dir + "JYD_20230112_BQ1.csv"
    # plink_file_230206_BQ2_fDR = plink_root_dir + "JYD_20230112_BQ2.csv"
    # plink_file_230206_BQ5_fDR = plink_root_dir + "JYD_20230112_BQ5.csv"
    # plink_file_230206_BQ10_fDR = plink_root_dir + "JYD_20230112_BQ10.csv"
    # plink_file_230206_BQ25_fDR = plink_root_dir + "JYD_20230112_BQ25.csv"
    # plink_file_230206_BQ50_fDR = plink_root_dir + "JYD_20230112_BQ50.csv"
    # plink_file_230206_BQ100_fDR = plink_root_dir + "JYD_20230112_BQ100.csv"
    # plink_file_230206_BQ200_fDR = plink_root_dir + "JYD_20230112_BQ200.csv"
    # plink_file_230206_BQ0_fDR = plink_root_dir + "JYD_20230112_BQ0.csv"

    # plink_root_dir = "G:/MSdata/221129BQM/BQ2/reports/crosslink_spotlinkform/"
    # plink_file_221121_BQ1_fDR = plink_root_dir + "JYD_20221121_BQ_1.csv"
    # plink_file_221121_BQ2_fDR = plink_root_dir + "JYD_20221121_BQ_2.csv"
    # plink_file_221121_BQ5_fDR = plink_root_dir + "JYD_20221121_BQ_5.csv"
    # plink_file_221121_BQ10_fDR = plink_root_dir + "JYD_20221121_BQ_10.csv"
    # plink_file_221121_BQ25_fDR = plink_root_dir + "JYD_20221121_BQ_25.csv"
    # plink_file_221121_BQ50_fDR = plink_root_dir + "JYD_20221121_BQ_50.csv"
    # plink_file_221121_BQ100_fDR = plink_root_dir + "JYD_20221121_BQ_100.csv"
    # plink_file_221121_BQ200_fDR = plink_root_dir + "JYD_20221121_BQ_200.csv"
    # plink_file_221121_BQ0_fDR = plink_root_dir + "JYD_20221121_BQ_300.csv"

    # BQ1_230206_p = plink_report_valid_link_sfDR(plink_file_230206_BQ1_fDR,
    #                                          bq_fasta)[0]
    # BQ2_230206_p = plink_report_valid_link_sfDR(plink_file_230206_BQ2_fDR,
    #                                          bq_fasta)[0]
    # BQ5_230206_p = plink_report_valid_link_sfDR(plink_file_230206_BQ5_fDR,
    #                                          bq_fasta)[0]
    # BQ10_230206_p = plink_report_valid_link_sfDR(plink_file_230206_BQ10_fDR,
    #                                          bq_fasta)[0]
    # BQ25_230206_p = plink_report_valid_link_sfDR(plink_file_230206_BQ25_fDR,
    #                                           bq_fasta)[0]
    # BQ50_230206_p = plink_report_valid_link_sfDR(plink_file_230206_BQ50_fDR,
    #                                           bq_fasta)[0]
    # BQ100_230206_p = plink_report_valid_link_sfDR(plink_file_230206_BQ100_fDR,
    #                                            bq_fasta)[0]
    # BQ200_230206_p = plink_report_valid_link_sfDR(plink_file_230206_BQ200_fDR,
    #                                            bq_fasta)[0]
    # BQ0_230206_p = plink_report_valid_link_sfDR(plink_file_230206_BQ0_fDR,
    #                                            bq_fasta)[0]

    # BQ1_221121_p = plink_report_valid_link_sfDR(plink_file_221121_BQ1_fDR,
    #                                          bq_fasta)[0]
    # BQ2_221121_p = plink_report_valid_link_sfDR(plink_file_221121_BQ2_fDR,
    #                                          bq_fasta)[0]
    # BQ5_221121_p = plink_report_valid_link_sfDR(plink_file_221121_BQ5_fDR,
    #                                          bq_fasta)[0]
    # BQ10_221121_p = plink_report_valid_link_sfDR(plink_file_221121_BQ10_fDR,
    #                                          bq_fasta)[0]
    # BQ25_221121_p = plink_report_valid_link_sfDR(plink_file_221121_BQ25_fDR,
    #                                           bq_fasta)[0]
    # BQ50_221121_p = plink_report_valid_link_sfDR(plink_file_221121_BQ50_fDR,
    #                                           bq_fasta)[0]
    # BQ100_221121_p = plink_report_valid_link_sfDR(plink_file_221121_BQ100_fDR,
    #                                            bq_fasta)[0]
    # BQ200_221121_p = plink_report_valid_link_sfDR(plink_file_221121_BQ200_fDR,
    #                                            bq_fasta)[0]
    # BQ0_221121_p = plink_report_valid_link_sfDR(plink_file_221121_BQ0_fDR,
    #                                            bq_fasta)[0]

    # # BQ1_repeat = cal_multi_list([BQ1_230206_p, BQ1_230206_p])
    # # BQ2_repeat = cal_multi_list([BQ2_230206_p, BQ2_230206_p])
    # # BQ5_repeat = cal_multi_list([BQ5_230206_p, BQ5_230206_p])
    # # BQ10_repeat = cal_multi_list([BQ10_230206_p, BQ10_230206_p])
    # # BQ25_repeat = cal_multi_list([BQ25_230206_p, BQ25_230206_p])
    # # BQ50_repeat = cal_multi_list([BQ50_230206_p, BQ50_230206_p])
    # # BQ100_repeat = cal_multi_list([BQ100_230206_p, BQ100_230206_p])
    # # BQ200_repeat = cal_multi_list([BQ200_230206_p, BQ200_230206_p])
    # # BQ0_repeat = cal_multi_list([BQ0_230206_p, BQ0_230206_p])

    # BQ1_repeat = cal_multi_list([BQ1_230206_p, BQ1_221121_p])
    # BQ2_repeat = cal_multi_list([BQ2_230206_p, BQ2_221121_p])
    # BQ5_repeat = cal_multi_list([BQ5_230206_p, BQ5_221121_p])
    # BQ10_repeat = cal_multi_list([BQ10_230206_p, BQ10_221121_p])
    # BQ25_repeat = cal_multi_list([BQ25_230206_p, BQ25_221121_p])
    # BQ50_repeat = cal_multi_list([BQ50_230206_p, BQ50_221121_p])
    # BQ100_repeat = cal_multi_list([BQ100_230206_p, BQ100_221121_p])
    # BQ200_repeat = cal_multi_list([BQ200_230206_p, BQ200_221121_p])
    # BQ0_repeat = cal_multi_list([BQ0_230206_p, BQ0_221121_p])

    # BQ_repeat_list = [BQ1_repeat, BQ2_repeat, BQ5_repeat, BQ10_repeat, BQ25_repeat, BQ50_repeat, BQ100_repeat,BQ200_repeat,BQ0_repeat]
    # print('|-A-|-C-|-D-|-E-|-F-|-G-|-H-|-I-|-K-|-L-|-M-|-N-|-P-|-Q-|-R-|-S-|-T-|-V-|-W-|-Y-|ALL|')
    # print('|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|---|')
    # for i in BQ_repeat_list:
    #     k = [
    #         j for j in i
    #         if bq_fasta[int(j[0]) -
    #                              1] == "K" or bq_fasta[int(j[1]) -
    #                                                             1] == "K"
    #     ]
    #     # cal_distance_pos_list(k,
    #     #                       fasta=bq_fasta,
    #     #                       pdb_file='G:/MSdata/pdb/conalbumin.pdb',
    #     #                       threshold=400,
    #     #                       print_all=False)
    #     report_animo_ratio(k,
    #                        pos_form=True,
    #                        fasta_file=bq_fasta,
    #                        tableform=True,
    #                        ratio_output=False,
    #                        threshold=250)
    # for i in BQ_repeat_list:
    #     k = [
    #         j for j in i
    #         if bq_fasta[int(j[0]) -
    #                              1] == "K" or bq_fasta[int(j[1]) -
    #                                                             1] == "K"
    #     ]
    #     cal_distance_pos_list(k,
    #                           fasta=bq_fasta,
    #                           pdb_file='G:/MSdata/pdb/conalbumin.pdb',
    #                           threshold=300,
    #                           print_all=False)
    #     # report_animo_ratio(k,
    #     #                    pos_form=True,
    #     #                    fasta_file=bq_fasta,
    #     #                    tableform=True,
    #     #                    ratio_output=False,
    #     #                    threshold=400)
    

    # print("*" * 40)
    # print("--Lactoferrin--" * 5)
    # print("*" * 40)
    # ## define part
    # plink_root_dir = "G:/MSdata/221129BQM/M/reports/crosslink_spotlinkform/"
    # ### fasta
    # lactoferrin_fasta = "CTISQPEWFKCRRWQWRMKKLGAPSITCVRRAFALECIRAIAEKKADAVTLDGGMVFEACRDPYKLRPVAAEIYGTKESPQTHYYAVAVVKKGSNFQLDQLQGRKSCHTGLGRSAGWIIPMGILRPYLSWTESLEPLQGAVAKFFSASCVPCIDRQAYPNLCQLCKGEGENQCACSSREPYFGYSGAFKCLQDGAGDVAFVKETTVFENLPEKADRDQYELLCLNNSRAPVDAFKECHLAQVPSHAVVARSVDGKEDLIWKLLSKAQEKFGKNKSRSFQLFGSPPGQRDLLFKDSALGFLRIPSKVDSALYLGSRYLTTLKNLRETAEEVKARYTRVVWCAVGPEEQKKCQQWSQQSGQNVTCATASTTDDCIVLVLKGEADALNLDGGYIYTAGKCGLVPVLAENRKSSKHSSLDCVLRPTEGYLAVAVVKKANEGLTWNSLKDKKSCHTAVDRTAGWNIPMGLIVNQTGSCAFDEFFSQSCAPGADPKSRLCALCAGDDQGLDKCVPNSKEKYYGYTGAFRCLAEDVGDVAFVKNDTVWENTNGESTADWAKNLNREDFRLLCLDGTRKPVTEAQSCHLAVAPNHAVVSRSDRAAHVKQVLLHQQALFGKNGKNCPDKFCLFKSETKNLLFNDNTECLAKLGGRPTYEEYLGTEYVTAIANLKKCSTSPLLEACAFLTR"
    # ### file_dir
    # plink_file_221121_LBM1_fDR = plink_root_dir + "JYD_20221121_LBM_1.csv"
    # plink_file_221121_LBM2_fDR = plink_root_dir + "JYD_20221121_LBM_2.csv"
    # plink_file_221121_LBM5_fDR = plink_root_dir + "JYD_20221121_LBM_5.csv"
    # plink_file_221121_LBM10_fDR = plink_root_dir + "JYD_20221121_LBM_10.csv"
    # plink_file_221121_LBM25_fDR = plink_root_dir + "JYD_20221121_LBM_25.csv"
    # plink_file_221121_LBM50_fDR = plink_root_dir + "JYD_20221121_LBM_50.csv"
    # plink_file_221121_LBM100_fDR = plink_root_dir + "JYD_20221121_LBM_100.csv"
    # plink_file_221121_LBM200_fDR = plink_root_dir + "JYD_20221121_LBM_200.csv"
    # plink_file_221121_LBM0_fDR = plink_root_dir + "JYD_20221121_LBM_300.csv"

    # plink_root_dir = "G:/MSdata/230303BQM/LBM/reports/crosslink_spotlinkform/"

    # plink_file_230112_LBM1_fDR = plink_root_dir + "JYD_20230112_LBM1.csv"
    # plink_file_230112_LBM2_fDR = plink_root_dir + "JYD_20230112_LBM2.csv"
    # plink_file_230112_LBM5_fDR = plink_root_dir + "JYD_20230112_LBM5.csv"
    # plink_file_230112_LBM10_fDR = plink_root_dir + "JYD_20230112_LBM10.csv"
    # plink_file_230112_LBM25_fDR = plink_root_dir + "JYD_20230112_LBM25.csv"
    # plink_file_230112_LBM50_fDR = plink_root_dir + "JYD_20230112_LBM50.csv"
    # plink_file_230112_LBM100_fDR = plink_root_dir + "JYD_20230112_LBM100.csv"
    # plink_file_230112_LBM200_fDR = plink_root_dir + "JYD_20230112_LBM200.csv"
    # plink_file_230112_LBM0_fDR = plink_root_dir + "JYD_20230112_LBM0.csv"

    # LBM0_221121_p = plink_report_valid_link_sfDR(plink_file_221121_LBM0_fDR,
    #                                          lactoferrin_fasta)[0]
    # LBM1_221121_p = plink_report_valid_link_sfDR(plink_file_221121_LBM1_fDR,
    #                                          lactoferrin_fasta)[0]
    # LBM2_221121_p = plink_report_valid_link_sfDR(plink_file_221121_LBM2_fDR,
    #                                          lactoferrin_fasta)[0]
    # LBM5_221121_p = plink_report_valid_link_sfDR(plink_file_221121_LBM5_fDR,
    #                                          lactoferrin_fasta)[0]
    # LBM10_221121_p = plink_report_valid_link_sfDR(plink_file_221121_LBM10_fDR,
    #                                          lactoferrin_fasta)[0]
    # LBM25_221121_p = plink_report_valid_link_sfDR(plink_file_221121_LBM25_fDR,
    #                                           lactoferrin_fasta)[0]
    # LBM50_221121_p = plink_report_valid_link_sfDR(plink_file_221121_LBM50_fDR,
    #                                           lactoferrin_fasta)[0]
    # LBM100_221121_p = plink_report_valid_link_sfDR(plink_file_221121_LBM100_fDR,
    #                                            lactoferrin_fasta)[0]
    # LBM200_221121_p = plink_report_valid_link_sfDR(plink_file_221121_LBM200_fDR,
    #                                            lactoferrin_fasta)[0]

    # LBM0_230112_p = plink_report_valid_link_sfDR(plink_file_230112_LBM0_fDR,
    #                                          lactoferrin_fasta)[0]
    # LBM1_230112_p = plink_report_valid_link_sfDR(plink_file_230112_LBM1_fDR,
    #                                          lactoferrin_fasta)[0]
    # LBM2_230112_p = plink_report_valid_link_sfDR(plink_file_230112_LBM2_fDR,
    #                                          lactoferrin_fasta)[0]
    # LBM5_230112_p = plink_report_valid_link_sfDR(plink_file_230112_LBM5_fDR,
    #                                          lactoferrin_fasta)[0]
    # LBM10_230112_p = plink_report_valid_link_sfDR(plink_file_230112_LBM10_fDR,
    #                                          lactoferrin_fasta)[0]
    # LBM25_230112_p = plink_report_valid_link_sfDR(plink_file_230112_LBM25_fDR,
    #                                           lactoferrin_fasta)[0]
    # LBM50_230112_p = plink_report_valid_link_sfDR(plink_file_230112_LBM50_fDR,
    #                                           lactoferrin_fasta)[0]
    # LBM100_230112_p = plink_report_valid_link_sfDR(plink_file_230112_LBM100_fDR,
    #                                            lactoferrin_fasta)[0]
    # LBM200_230112_p = plink_report_valid_link_sfDR(plink_file_230112_LBM200_fDR,
    #                                            lactoferrin_fasta)[0]

    # LBM1_repeat = cal_multi_list([LBM1_230112_p, LBM1_221121_p])
    # LBM2_repeat = cal_multi_list([LBM2_230112_p, LBM2_221121_p])
    # LBM5_repeat = cal_multi_list([LBM5_230112_p, LBM5_221121_p])
    # LBM10_repeat = cal_multi_list([LBM10_230112_p, LBM10_221121_p])
    # LBM25_repeat = cal_multi_list([LBM25_230112_p, LBM25_221121_p])
    # LBM50_repeat = cal_multi_list([LBM50_230112_p, LBM50_221121_p])
    # LBM100_repeat = cal_multi_list([LBM100_230112_p, LBM100_221121_p])
    # LBM200_repeat = cal_multi_list([LBM200_230112_p, LBM200_221121_p])
    # LBM0_repeat = cal_multi_list([LBM0_230112_p, LBM0_221121_p])

    # # LBM1_repeat = cal_multi_list([LBM1_230112_p, LBM1_230112_p])
    # # LBM2_repeat = cal_multi_list([LBM2_230112_p, LBM2_230112_p])
    # # LBM5_repeat = cal_multi_list([LBM5_230112_p, LBM5_230112_p])
    # # LBM10_repeat = cal_multi_list([LBM10_230112_p, LBM10_230112_p])
    # # LBM25_repeat = cal_multi_list([LBM25_230112_p, LBM25_230112_p])
    # # LBM50_repeat = cal_multi_list([LBM50_230112_p, LBM50_230112_p])
    # # LBM100_repeat = cal_multi_list([LBM100_230112_p, LBM100_230112_p])
    # # LBM200_repeat = cal_multi_list([LBM200_230112_p, LBM200_230112_p])
    # # LBM0_repeat = cal_multi_list([LBM0_230112_p, LBM0_230112_p])


    # M_repeat_list = [LBM0_repeat, LBM1_repeat, LBM2_repeat, LBM5_repeat, LBM10_repeat, LBM25_repeat, LBM50_repeat,LBM100_repeat,LBM200_repeat]
    # print('|-A-|-C-|-D-|-E-|-F-|-G-|-H-|-I-|-K-|-L-|-M-|-N-|-P-|-Q-|-R-|-S-|-T-|-V-|-W-|-Y-|ALL|')
    # print('|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|---|')
    # for i in M_repeat_list:
    #     k = [
    #         j for j in i
    #         if lactoferrin_fasta[int(j[0]) -
    #                              1] == "K" or lactoferrin_fasta[int(j[1]) -
    #                                                             1] == "K"
    #     ]
    #     # cal_distance_pos_list(k,
    #     #                       fasta=lactoferrin_fasta,
    #     #                       pdb_file='G:/MSdata/pdb/lactoferrin.pdb',
    #     #                       threshold=40,
    #     #                       print_all=False)
    #     report_animo_ratio(k,
    #                        pos_form=True,
    #                        fasta_file=lactoferrin_fasta,
    #                        tableform=True,
    #                        ratio_output=False,
    #                        threshold=40)
        
    # for i in M_repeat_list:
    #     k = [
    #         j for j in i
    #         if lactoferrin_fasta[int(j[0]) -
    #                              1] == "K" or lactoferrin_fasta[int(j[1]) -
    #                                                             1] == "K"
    #     ]
    #     cal_distance_pos_list(k,
    #                           fasta=lactoferrin_fasta,
    #                           pdb_file='G:/MSdata/pdb/lactoferrin.pdb',
    #                           threshold=40,
    #                           print_all=False)
    #     # report_animo_ratio(k,
    #     #                    pos_form=True,
    #     #                    fasta_file=lactoferrin_fasta,
    #     #                    tableform=True,
    #     #                    ratio_output=False,
    #     #                    threshold=40)

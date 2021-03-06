import csv
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
            print("Find pep Fail:", pep)
            return 0, 0


def find_link_pos(ms_peptide, fasta):
    pep1, pep2 = ms_peptide.split('-')
    pep1_fasta, pep1_num = pep1.strip(')').split('(')
    pep2_fasta, pep2_num = pep2.strip(')').split('(')
    cross_pos_1 = pep1_fasta[int(pep1_num) - 1]
    cross_pos_2 = pep2_fasta[int(pep2_num) - 1]
    pep1_start, pep1_end = find_loc_in_fasta(pep1_fasta, fasta)
    pep2_start, pep2_end = find_loc_in_fasta(pep2_fasta,
                                             fasta)  # if not find return 0,0

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


### Analyze amino ratio
# Input list like [["K","L"],["Y","K"],??????]
# if pos_form on, input should be [[139,145],[224,269],??????]


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
    else:
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

    # Report
    if tableform == False:
        return cross_dic
    # Out dic like {'A': 31, 'C': 26, 'D': 72, 'E': 117, 'F': 11, 'G': 25,??????
    else:
        print('|A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y|ALL|')
        print('|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|---|')
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
    # Out dic like {'A': 31, 'C': 26, 'D': 72, 'E': 117, 'F': 11, 'G': 25,??????
    else:
        print('|A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y|ALL|')
        print('|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|---|')
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
        least_ssn=2,
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
            if int(i[2][2]) >= least_ssn
        ]
    elif filtermode == True and f2sf == False:
        report_list = [
            i for i in sorted(report_list, key=lambda x: x[2][1], reverse=True)
            [0:tsd]
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
        least_ssn=1,
        threshold=0.8,
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
    atom1 = residue1[animo_core1]
    atom2 = residue2[animo_core2]
    distance = atom1 - atom2
    if autoprint == True:
        print('Distance between', residue1.get_resname(), animo_core1, 'and',
              residue2.get_resname(), animo_core2, 'is ', distance)
    return distance


def cal_distance_pos_list(
        pos_list,
        fasta="MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA",
        pdb_file='G:/MSdata/bsa.pdb',
        max_distance=20,
        min_distance=6,
        read_mod=False):

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
        print(i, dis, fasta[int(i[0]) - 1], fasta[int(i[1]) - 1])
        dis_list[0].append(i)
        dis_list[1].append(dis)

    print("Average", format(np.mean(dis_list[1]), '.5f'), ', Median',
          format(np.median(dis_list[1]), '.5f'), ', Std', np.std(dis_list[1]),
          ', Num', len(dis_list[1]), ', Over Max Num', over_max_num,
          ', Over Min Num', over_min_num)
    if read_mod == True:
        dis_list = list(map(list, zip(*dis_list)))

        return dis_list
    else:
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
        print('|A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y|ALL|')
        print('|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|---|')
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
                        cross_dic["all"] += 1
                except:
                    continue
    if table_form == True:
        print('|A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y|ALL|')
        print('|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|---|')
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
    atom1 = residue1[animo_core1]
    atom2 = residue2[animo_core2]
    distance = atom1 - atom2
    if autoprint == True:
        print('Distance between', residue1.get_resname(), animo_core1, 'and',
              residue2.get_resname(), animo_core2, 'is ', distance)
    return distance


# if __name__ == '__main__':

bsa_fasta = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA"

# plink2crosslink(
#     "G:/MSdata/220321BSA_modify_power/PLINK/pLink_task_2022.03.26.09.56.22/reports/bsa_2022.03.26.filtered_cross-linked_spectra.csv",
#     "G:/MSdata/220321BSA_modify_power/PLINK/pLink_task_2022.03.26.09.56.22/reports/crosslink_spotlinkform/"
# )

spotlink_file_0322_B_fDR = 'G:/MSdata/220321BSA_modify_power/SPOTLINK/JYD_20220322_B_HCDFT_result_filtered.csv'
spotlink_file_0322_BO_fDR = 'G:/MSdata/220321BSA_modify_power/SPOTLINK/JYD_20220322_BO_HCDFT_result_filtered.csv'
spotlink_file_0322_BN_fDR = 'G:/MSdata/220321BSA_modify_power/SPOTLINK/JYD_20220322_BN_HCDFT_result_filtered.csv'
spotlink_file_0322_BU_fDR = 'G:/MSdata/220321BSA_modify_power/SPOTLINK/JYD_20220322_BU_HCDFT_result_filtered.csv'
spotlink_file_0322_B1_fDR = 'G:/MSdata/220321BSA_modify_power/SPOTLINK/JYD_20220322_B1_HCDFT_result_filtered.csv'
spotlink_file_0322_B5_fDR = 'G:/MSdata/220321BSA_modify_power/SPOTLINK/JYD_20220322_B5_HCDFT_result_filtered.csv'
spotlink_file_0322_B50_fDR = 'G:/MSdata/220321BSA_modify_power/SPOTLINK/JYD_20220322_B50_HCDFT_result_filtered.csv'
spotlink_file_0322_B500_fDR = 'G:/MSdata/220321BSA_modify_power/SPOTLINK/JYD_20220322_B500_HCDFT_result_filtered.csv'
spotlink_file_0322_B1000_fDR = 'G:/MSdata/220321BSA_modify_power/SPOTLINK/JYD_20220322_B1000_HCDFT_result_filtered.csv'

crosslink_file_0118_BM10_fDR = 'G:/MSdata/220118UBBSA/spotlink/WXZ_20220118_BM10_HCDFT_result_filtered.csv'

plink_0322_B_fDR = "G:/MSdata/220321BSA_modify_power/PLINK/pLink_task_2022.03.26.09.56.22/reports/crosslink_spotlinkform/JYD_20220322_B.csv"
plink_0322_BO_fDR = "G:/MSdata/220321BSA_modify_power/PLINK/pLink_task_2022.03.26.09.56.22/reports/crosslink_spotlinkform/JYD_20220322_BO.csv"
plink_0322_BU_fDR = "G:/MSdata/220321BSA_modify_power/PLINK/pLink_task_2022.03.26.09.56.22/reports/crosslink_spotlinkform/JYD_20220322_BU.csv"
plink_0322_BN_fDR = "G:/MSdata/220321BSA_modify_power/PLINK/pLink_task_2022.03.26.09.56.22/reports/crosslink_spotlinkform/JYD_20220322_BN.csv"
plink_0322_B1_fDR = "G:/MSdata/220321BSA_modify_power/PLINK/pLink_task_2022.03.26.09.56.22/reports/crosslink_spotlinkform/JYD_20220322_B1.csv"
plink_0322_B5_fDR = "G:/MSdata/220321BSA_modify_power/PLINK/pLink_task_2022.03.26.09.56.22/reports/crosslink_spotlinkform/JYD_20220322_B5.csv"
plink_0322_B50_fDR = "G:/MSdata/220321BSA_modify_power/PLINK/pLink_task_2022.03.26.09.56.22/reports/crosslink_spotlinkform/JYD_20220322_B50.csv"
plink_0322_B500_fDR = "G:/MSdata/220321BSA_modify_power/PLINK/pLink_task_2022.03.26.09.56.22/reports/crosslink_spotlinkform/JYD_20220322_B500.csv"
plink_0322_B1000_fDR = "G:/MSdata/220321BSA_modify_power/PLINK/pLink_task_2022.03.26.09.56.22/reports/crosslink_spotlinkform/JYD_20220322_B1000.csv"

B_0322_s = spotlink_report_valid_link_sfDR(spotlink_file_0322_B_fDR,
                                           bsa_fasta)[0]
BO_0322_s = spotlink_report_valid_link_sfDR(spotlink_file_0322_BO_fDR,
                                            bsa_fasta)[0]
BN_0322_s = spotlink_report_valid_link_sfDR(spotlink_file_0322_BN_fDR,
                                            bsa_fasta)[0]
BU_0322_s = spotlink_report_valid_link_sfDR(spotlink_file_0322_BU_fDR,
                                            bsa_fasta)[0]
B1_0322_s = spotlink_report_valid_link_sfDR(spotlink_file_0322_B1_fDR,
                                            bsa_fasta)[0]
B5_0322_s = spotlink_report_valid_link_sfDR(spotlink_file_0322_B5_fDR,
                                            bsa_fasta)[0]
B50_0322_s = spotlink_report_valid_link_sfDR(spotlink_file_0322_B50_fDR,
                                             bsa_fasta)[0]
B500_0322_s = spotlink_report_valid_link_sfDR(spotlink_file_0322_B500_fDR,
                                              bsa_fasta)[0]
B1000_0322_s = spotlink_report_valid_link_sfDR(spotlink_file_0322_B1000_fDR,
                                               bsa_fasta)[0]
BM10_0118_s = spotlink_report_valid_link_sfDR(crosslink_file_0118_BM10_fDR,
                                              bsa_fasta,
                                              least_ssn=1,
                                              threshold=0.8)[0]

B_0322_p = plink_report_valid_link_sfDR(plink_0322_B_fDR, bsa_fasta)[0]
BO_0322_p = plink_report_valid_link_sfDR(plink_0322_BO_fDR, bsa_fasta)[0]
BU_0322_p = plink_report_valid_link_sfDR(plink_0322_BU_fDR, bsa_fasta)[0]
BN_0322_p = plink_report_valid_link_sfDR(plink_0322_BN_fDR, bsa_fasta)[0]
B1_0322_p = plink_report_valid_link_sfDR(plink_0322_B1_fDR, bsa_fasta)[0]
B5_0322_p = plink_report_valid_link_sfDR(plink_0322_B5_fDR, bsa_fasta)[0]
B50_0322_p = plink_report_valid_link_sfDR(plink_0322_B50_fDR, bsa_fasta)[0]
B500_0322_p = plink_report_valid_link_sfDR(plink_0322_B500_fDR, bsa_fasta)[0]
B1000_0322_p = plink_report_valid_link_sfDR(plink_0322_B1000_fDR, bsa_fasta)[0]
print("-------------")

B_repeat = cal_repeat_list(B_0322_p, B_0322_s)[0]
BO_repeat = cal_repeat_list(BO_0322_p, BO_0322_s)[0]
BU_repeat = cal_repeat_list(BU_0322_p, BU_0322_s)[0]
BN_repeat = cal_repeat_list(BN_0322_p, BN_0322_s)[0]
B1_repeat = cal_repeat_list(B1_0322_p, B1_0322_s)[0]
B5_repeat = cal_repeat_list(B5_0322_p, B5_0322_s)[0]
B50_repeat = cal_repeat_list(B50_0322_p, B50_0322_s)[0]
B500_repeat = cal_repeat_list(B500_0322_p, B500_0322_s)[0]
B1000_repeat = cal_repeat_list(B1000_0322_p, B1000_0322_s)[0]
# Contrast with standard
print("-------------")
print(len(BM10_0118_s))
B_repeat = cal_repeat_list(B_repeat, BM10_0118_s)[0]
B1_repeat = cal_repeat_list(B1_repeat, BM10_0118_s)[0]
BO_repeat = cal_repeat_list(BO_repeat, BM10_0118_s)[0]
B5_repeat = cal_repeat_list(B5_repeat, BM10_0118_s)[0]
B50_repeat = cal_repeat_list(B50_repeat, BM10_0118_s)[0]
B500_repeat = cal_repeat_list(B500_repeat, BM10_0118_s)[0]
B1000_repeat = cal_repeat_list(B1000_repeat, BM10_0118_s)[0]
print("-------------")

repeat_list_0322 = [
    B_repeat, B1_repeat, BO_repeat, B5_repeat, B50_repeat, B500_repeat,
    B1000_repeat
]

list_0322_s = [
    B_0322_s, BO_0322_s, BN_0322_s, BU_0322_s, B1_0322_s, B5_0322_s,
    B50_0322_s, B500_0322_s, B1000_0322_s
]

list_0322_p = [
    B_0322_p, BO_0322_p, BN_0322_p, BU_0322_p, B1_0322_p, B5_0322_p,
    B50_0322_p, B500_0322_p, B1000_0322_p
]

for i in repeat_list_0322:
    report_animo_ratio(i,
                       pos_form=True,
                       fasta_file=bsa_fasta,
                       tableform=True,
                       ratio_output=True,
                       threshold=40)

# adk_fasta = "MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG"
# # plink2crosslink(
# #     "G:/MSdata/220329ADK/adk/plink/pLink_task_2022.04.02.20.19.24/reports/adk_2022.04.02.filtered_cross-linked_spectra.csv",
# #     "G:/MSdata/220329ADK/adk/plink/pLink_task_2022.04.02.20.19.24/reports/crosslink_spotlinkform/"
# # )
# plink_0329_ADK_fDR = "G:/MSdata/220329ADK/adk/plink/pLink_task_2022.04.02.20.19.24/reports/crosslink_spotlinkform/WXZ_20220329_AdK.csv"
# plink_0329_ADKC_fDR = "G:/MSdata/220329ADK/adk/plink/pLink_task_2022.04.02.20.19.24/reports/crosslink_spotlinkform/WXZ_20220329_AdKC.csv"

# spotlink_file_0329_ADK_fDR = 'G:/MSdata/220329ADK/adk/spotlink/WXZ_20220329_AdK_HCDFT_result_filtered.csv'
# spotlink_file_0329_ADKC_fDR = 'G:/MSdata/220329ADK/adk/spotlink/WXZ_20220329_AdKC_HCDFT_result_filtered.csv'

# ADK_0329_s = spotlink_report_valid_link_sfDR(spotlink_file_0329_ADK_fDR,
#                                              adk_fasta)[0]
# ADKC_0329_s = spotlink_report_valid_link_sfDR(spotlink_file_0329_ADKC_fDR,
#                                               adk_fasta)[0]

# ADK_0329_p = plink_report_valid_link_sfDR(plink_0329_ADK_fDR, adk_fasta)[0]
# ADKC_0329_p = plink_report_valid_link_sfDR(plink_0329_ADKC_fDR, adk_fasta)[0]

# ADK_repeat = cal_repeat_list(ADK_0329_s, ADK_0329_p)[0]
# ADKC_repeat = cal_repeat_list(ADKC_0329_s, ADKC_0329_p)[0]

# report_animo_ratio(ADK_repeat,
#                    pos_form=True,
#                    fasta_file=adk_fasta,
#                    tableform=True)
# report_animo_ratio(ADKC_repeat,
#                    pos_form=True,
#                    fasta_file=adk_fasta,
#                    tableform=True)
# calc_fasta(adk_fasta, table_form=True)

# from pymol import cmd

# for i in ADKC_repeat:
#     dis = cal_distance(int(i[0]), int(i[1]), 'G:/MSdata/220329ADK/adk/adk.pdb')
#     print(i, dis, adk_fasta[int(i[0]) - 1], adk_fasta[int(i[1]) - 1])
#     cmd.distance(
#         str(i[0]) + '-' + str(i[1]), "(/adk//A/" + str(i[0]) + "/CA)",
#         "( /adk//A/" + str(i[1]) + "/CA)")

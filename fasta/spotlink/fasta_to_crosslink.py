import csv
import pdb_distance_analyze
import numpy as np


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


def report_link_pos(ms_csv_file, fasta, ms_peptide_list=5):
    # Spotlink
    link_pos_list = []
    with open(ms_csv_file, "r") as f:
        reader = csv.reader(f)
        for target in reader:
            try:
                link_pos_list.append(
                    find_link_pos(target[ms_peptide_list], fasta))

                link_pos_list[-1].append([
                    target[ms_peptide_list + 5], target[ms_peptide_list + 6],
                    target[ms_peptide_list + 7]
                ])
            except:
                if target[
                        ms_peptide_list] != "Alpha Peptide Protein" and target[
                            ms_peptide_list] != "Peptide":
                    print("Target link position cannot be analyzed:",
                          target[ms_peptide_list])
    return link_pos_list


# this function will del repeat element
def cal_repeat_list(list_1, list_2, simple_mode=False):
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


def fdr2sfdr(file_name):
    # Spotlink
    return 'site'.join(file_name.split('result'))


def cal_same_link_pos(ms_file_1,
                      ms_file_2,
                      fasta,
                      set=False,
                      simple_mode=False,
                      filtermode=True,
                      threshold=0.75,
                      ms_peptide_list=5,
                      f2sf=False):
    #Spotlink
    if f2sf == True:
        ms_file_1 = fdr2sfdr(ms_file_1)
        ms_file_2 = fdr2sfdr(ms_file_2)
        ms_peptide_list = 2
    report_list_1 = report_link_pos(ms_file_1,
                                    fasta,
                                    ms_peptide_list=ms_peptide_list)
    report_list_2 = report_link_pos(ms_file_2,
                                    fasta,
                                    ms_peptide_list=ms_peptide_list)

    if filtermode == True:
        report_list_1 = sorted(report_list_1,
                               key=lambda x: x[2][1],
                               reverse=True)[0:int(threshold *
                                                   len(report_list_1))]
        report_list_2 = sorted(report_list_2,
                               key=lambda x: x[2][1],
                               reverse=True)[0:int(threshold *
                                                   len(report_list_2))]

    if set == False:
        hard_loc_list_1 = [[i[0][-1], i[1][-1]] for i in report_list_1]
        hard_loc_list_2 = [[i[0][-1], i[1][-1]] for i in report_list_2]
    elif set == True:
        hard_loc_list_1 = [[
            max([i[0][-1], i[1][-1]]),
            min([i[0][-1], i[1][-1]])
        ] for i in report_list_1]
        hard_loc_list_2 = [[
            max([i[0][-1], i[1][-1]]),
            min([i[0][-1], i[1][-1]])
        ] for i in report_list_2]
    repeat_list, unrepeat_list_1, unrepeat_list_2, repeat_ratio_1, repeat_ratio_2 = cal_repeat_list(
        hard_loc_list_1, hard_loc_list_2)
    if simple_mode == False:
        return [
            repeat_list, unrepeat_list_1, unrepeat_list_2, repeat_ratio_1,
            repeat_ratio_2
        ]
    else:
        print(repeat_ratio_1, repeat_ratio_2, len(repeat_list),
              len(hard_loc_list_1), len(hard_loc_list_2))
        report_animo_ratio(repeat_list,
                           pos_form=True,
                           fasta_file=fasta,
                           tableform=True)
        return [repeat_list]


def report_animo_ratio(
        link_list,
        hard_animo="K",
        pos_form=False,
        fasta_file='MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA',
        tableform=False):
    # Input list like [["K","L"],["Y","K"],??????]
    if pos_form == True:
        pos_list = link_list
        link_list = []
        for i in pos_list:
            link_list.append(
                [fasta_file[int(i[0]) - 1], fasta_file[int(i[1]) - 1]])
    animo_list = list("ACDEFGHIKLMNPQRSTVWY")
    cross_dic = {}
    for animo_tag in animo_list:
        cross_dic[animo_tag] = 0
    cross_dic["all"] = 0
    cross_dic["nofind"] = 0

    for cross_dimer in link_list:
        if cross_dimer[0] == hard_animo:
            cross_dic[cross_dimer[1]] += 1
            cross_dic["all"] += 1
        elif cross_dimer[1] == hard_animo:
            cross_dic[cross_dimer[0]] += 1
            cross_dic["all"] += 1
        else:
            cross_dic["nofind"] += 0
    if tableform == False:
        return cross_dic
    # Out dic like {'A': 31, 'C': 26, 'D': 72, 'E': 117, 'F': 11, 'G': 25,??????
    else:
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


# Open spotlink file and report ms_peptide
def draw_cross_pep_list(dir, ms_peptide_list=5):
    cross_pep_list = []
    with open(dir) as f:
        reader = csv.reader(f)
        for target in reader:
            if "-" in target[ms_peptide_list]:
                cross_pep_list.append(target[ms_peptide_list])
    # print(cross_pep_list)
    return cross_pep_list


def report_animo_dic(ms_csv_dir, hard_animo="K"):
    cross_list = [
        [i[0][2], i[1][2]] for i in
        [find_link_pos(j, bsa_fasta) for j in draw_cross_pep_list(ms_csv_dir)]
        if 0 not in [i[0][3], i[0][4], i[1][3], i[1][4]]
    ]
    return report_animo_ratio(cross_list, hard_animo=hard_animo)


def count_element(target_list):
    count_list = [[], []]
    for i in target_list:
        if i not in count_list[0]:
            count_list[0].append(i)
            count_list[1].append(target_list.count(i))

    return list(map(list, zip(*count_list)))


def filter_baseonrepeat(target_list, repeat_times=3):
    repeat_list = count_element(target_list)
    print(repeat_list)
    filtered_list = []

    for i in repeat_list:
        if i[1] >= repeat_times:
            filtered_list.append(i[0])

    print('All elements', len(repeat_list), ', repeat threshold', repeat_times,
          ', Reserved elements', len(filtered_list))

    return filtered_list


def cal_distance(pos_list,
                 pdb_file='G:/MSdata/bsa.pdb',
                 max_distance=16,
                 min_distance=5):

    dis_list = [[], []]
    over_max_num = 0
    over_min_num = 0
    for i in pos_list:
        dis = float(
            format(
                pdb_distance_analyze.cal_distance(int(i[0]), int(i[1]),
                                                  pdb_file), '.5f'))
        if dis >= max_distance:
            over_max_num += 1
        elif dis <= min_distance:
            over_min_num += 1
        print(i, dis, bsa_fasta[int(i[0]) - 1], bsa_fasta[int(i[1]) - 1])
        dis_list[0].append(i)
        dis_list[1].append(dis)

    print("Average", format(np.mean(dis_list[1]), '.5f'), ', Median',
          format(np.median(dis_list[1]), '.5f'), ', Num', len(dis_list[1]),
          ', Over Max Num', over_max_num, ', Over Min Num', over_min_num)

    dis_list = list(map(list, zip(*dis_list)))

    return dis_list


def report_valid_link_sfDR(
        ms_csv_dir,
        fasta="MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA",
        set=True,
        filtermode=True,
        least_ssn=2,
        threshold=0.75,
        ms_peptide_list=2,
        f2sf=True):
    # Spotlink
    if f2sf == True:
        ms_csv_dir = fdr2sfdr(ms_csv_dir)
        ms_peptide_list = 2

    report_list = report_link_pos(ms_csv_dir,
                                  fasta,
                                  ms_peptide_list=ms_peptide_list)

    if filtermode == True:
        report_list = [
            i for i in sorted(report_list, key=lambda x: x[2][1], reverse=True)
            [0:int(threshold * len(report_list))] if int(i[2][2]) >= least_ssn
        ]

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

    report_animo_ratio(hard_loc_list,
                       pos_form=True,
                       fasta_file=fasta,
                       tableform=True)
    return [hard_loc_list]


if __name__ == "__main__":
    # Def
    bsa_fasta = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA"

    crosslink_file_1108_M10_fDR = 'G:/MSdata/211200BSA/20220120-plink/WXZ_20211108_60_10_HCDFT_result_filtered.csv'
    crosslink_file_1108_M30_fDR = 'G:/MSdata/211200BSA/20220120-plink/WXZ_20211108_60_30_HCDFT_result_filtered.csv'
    crosslink_file_1108_M15_fDR = 'G:/MSdata/211200BSA/20220120-plink/WXZ_20211108_60_15_HCDFT_result_filtered.csv'
    crosslink_file_1108_M3_fDR = 'G:/MSdata/211200BSA/20220120-plink/WXZ_20211108_60_3_HCDFT_result_filtered.csv'
    crosslink_file_1108_N_fDR = 'G:/MSdata/211200BSA/20220120-plink/WXZ_20211108_BSA_HCDFT_result_filtered.csv'

    crosslink_file_0114_M10_fDR = 'G:/MSdata/220118BSA/20220120-plink/JYD_20220114_M10_HCDFT_result_filtered.csv'
    crosslink_file_0114_M30_fDR = 'G:/MSdata/220118BSA/20220120-plink/JYD_20220114_M30_HCDFT_result_filtered.csv'
    crosslink_file_0114_M1_fDR = 'G:/MSdata/220118BSA/20220120-plink/JYD_20220114_M1_HCDFT_result_filtered.csv'
    crosslink_file_0114_S30_fDR = 'G:/MSdata/220118BSA/20220120-plink/JYD_20220114_S30_HCDFT_result_filtered.csv'
    crosslink_file_0114_S15_fDR = 'G:/MSdata/220118BSA/20220120-plink/JYD_20220114_S15_HCDFT_result_filtered.csv'
    crosslink_file_0114_S10_fDR = 'G:/MSdata/220118BSA/20220120-plink/JYD_20220114_S10_HCDFT_result_filtered.csv'
    crosslink_file_0114_S5_fDR = 'G:/MSdata/220118BSA/20220120-plink/JYD_20220114_S5_HCDFT_result_filtered.csv'
    crosslink_file_0114_S1_fDR = 'G:/MSdata/220118BSA/20220120-plink/JYD_20220114_S1_HCDFT_result_filtered.csv'
    crosslink_file_0114_N_fDR = 'G:/MSdata/220118BSA/20220120-plink/JYD_20220114_N_HCDFT_result_filtered.csv'

    crosslink_file_0118_B_fDR = 'G:/MSdata/220125UBBSA/spotlink/WXZ_20220118_B_HCDFT_result_filtered.csv'
    crosslink_file_0118_BM10_fDR = 'G:/MSdata/220125UBBSA/spotlink/WXZ_20220118_BM10_HCDFT_result_filtered.csv'
    crosslink_file_0118_BS1_fDR = 'G:/MSdata/220125UBBSA/spotlink/WXZ_20220118_BS1_HCDFT_result_filtered.csv'
    crosslink_file_0118_BS2_fDR = 'G:/MSdata/220125UBBSA/spotlink/WXZ_20220118_BS2_HCDFT_result_filtered.csv'
    crosslink_file_0118_BS3_fDR = 'G:/MSdata/220125UBBSA/spotlink/WXZ_20220118_BS3_HCDFT_result_filtered.csv'
    crosslink_file_0118_BS5_fDR = 'G:/MSdata/220125UBBSA/spotlink/WXZ_20220118_BS5_HCDFT_result_filtered.csv'
    crosslink_file_0118_BS10_fDR = 'G:/MSdata/220125UBBSA/spotlink/WXZ_20220118_BS10_HCDFT_result_filtered.csv'
    crosslink_file_0118_BS30_fDR = 'G:/MSdata/220125UBBSA/spotlink/WXZ_20220118_BS30_HCDFT_result_filtered.csv'

    BM10_0118 = report_valid_link_sfDR(crosslink_file_0118_BM10_fDR,
                                       bsa_fasta)[0]
    BM10_0114 = report_valid_link_sfDR(crosslink_file_0114_M10_fDR,
                                       bsa_fasta)[0]
    BM10_1108 = report_valid_link_sfDR(crosslink_file_1108_M10_fDR,
                                       bsa_fasta)[0]
    BM15_1108 = report_valid_link_sfDR(crosslink_file_1108_M15_fDR,
                                       bsa_fasta)[0]
    BM10_all = BM10_0118 + BM10_0114
    fBM10_all = filter_baseonrepeat(BM10_all, repeat_times=2)

    BS5_0118 = report_valid_link_sfDR(crosslink_file_0118_BS5_fDR,
                                      bsa_fasta)[0]
    BS5_0114 = report_valid_link_sfDR(crosslink_file_0114_S5_fDR, bsa_fasta)[0]
    BS10_0118 = report_valid_link_sfDR(crosslink_file_0118_BS10_fDR,
                                       bsa_fasta)[0]
    BS10_0114 = report_valid_link_sfDR(crosslink_file_0114_S10_fDR,
                                       bsa_fasta)[0]
    BS3_0118 = report_valid_link_sfDR(crosslink_file_0118_BS3_fDR,
                                      bsa_fasta)[0]
    BS2_0118 = report_valid_link_sfDR(crosslink_file_0118_BS2_fDR,
                                      bsa_fasta)[0]
    BS1_0118 = report_valid_link_sfDR(crosslink_file_0118_BS1_fDR,
                                      bsa_fasta)[0]

    BS5_all = BS5_0118 + BS5_0114

    fBS5_all = filter_baseonrepeat(BS5_all, repeat_times=2)
    cal_distance(fBM10_all)

    # repeat_list_S30 = cal_same_link_pos(crosslink_file_0118_BS30_fDR,
    #                                     crosslink_file_0118_BM10_fDR,
    #                                     bsa_fasta,
    #                                     set=True,
    #                                     simple_mode=True,
    #                                     f2sf=True)[0]

    # repeat_list_S10 = cal_same_link_pos(crosslink_file_0118_BS10_fDR,
    #                                     crosslink_file_0118_BM10_fDR,
    #                                     bsa_fasta,
    #                                     set=True,
    #                                     simple_mode=True,
    #                                     f2sf=False)[0]

    # repeat_list_S5 = cal_same_link_pos(crosslink_file_0118_BS5_fDR,
    #                                    crosslink_file_0118_BM10_fDR,
    #                                    bsa_fasta,
    #                                    set=True,
    #                                    simple_mode=True,
    #                                    f2sf=True)[0]

    # repeat_list_S3 = cal_same_link_pos(crosslink_file_0118_BS3_fDR,
    #                                    crosslink_file_0118_BM10_fDR,
    #                                    bsa_fasta,
    #                                    set=True,
    #                                    simple_mode=True,
    #                                    f2sf=True)[0]

    # repeat_list_S2 = cal_same_link_pos(crosslink_file_0118_BS2_fDR,
    #                                    crosslink_file_0118_BM10_fDR,
    #                                    bsa_fasta,
    #                                    set=True,
    #                                    simple_mode=True,
    #                                    f2sf=True)[0]

    # repeat_list_S1 = cal_same_link_pos(crosslink_file_0118_BS1_fDR,
    #                                    crosslink_file_0118_BM10_fDR,
    #                                    bsa_fasta,
    #                                    set=True,
    #                                    simple_mode=True,
    #                                    f2sf=True)[0]

    # repeat_list_N = cal_same_link_pos(crosslink_file_0118_B_fDR,
    #                                   crosslink_file_0118_BM10_fDR,
    #                                   bsa_fasta,
    #                                   set=True,
    #                                   simple_mode=True,
    #                                   f2sf=True)[0]

    # a = cal_repeat_list(repeat_list_S30, repeat_list_S10, simple_mode=True)[0]
    # a = cal_repeat_list(repeat_list_S30, repeat_list_S5, simple_mode=True)[0]
    # a = cal_repeat_list(repeat_list_S30, repeat_list_S3, simple_mode=True)[0]
    # a = cal_repeat_list(repeat_list_S30, repeat_list_S2, simple_mode=True)[0]
    # a = cal_repeat_list(repeat_list_S30, repeat_list_S1, simple_mode=True)[0]
    # a = cal_repeat_list(repeat_list_S30, repeat_list_N, simple_mode=True)[0]

    # import pdb_distance_analyze

    # for i in a:
    #     print(
    #         i,
    #         pdb_distance_analyze.cal_distance(int(i[0]), int(i[1]),
    #                                           'G:/MSdata/bsa.pdb'),
    #         bsa_fasta[int(i[0]) - 1], bsa_fasta[int(i[1]) - 1])

import csv
from os import uname_result
# Def
bsa_fasta = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA"

crosslink_file_30_10 = 'C:/Users/jiang/OneDrive/科研/tc/实验记录/2021/20211105 365crosslink/20211115_tangchun_BSA/WXZ_20211108_30_10_HCDFT_result_filtered.csv'
crosslink_file_60_3 = 'C:/Users/jiang/OneDrive/科研/tc/实验记录/2021/20211105 365crosslink/20211115_tangchun_BSA/WXZ_20211108_60_3_HCDFT_result_filtered.csv'
crosslink_file_60_10 = 'C:/Users/jiang/OneDrive/科研/tc/实验记录/2021/20211105 365crosslink/20211115_tangchun_BSA/WXZ_20211108_60_10_HCDFT_result_filtered.csv'
crosslink_file_60_15 = 'C:/Users/jiang/OneDrive/科研/tc/实验记录/2021/20211105 365crosslink/20211115_tangchun_BSA/WXZ_20211108_60_15_HCDFT_result_filtered.csv'
crosslink_file_60_30 = 'C:/Users/jiang/OneDrive/科研/tc/实验记录/2021/20211105 365crosslink/20211115_tangchun_BSA/WXZ_20211108_60_30_HCDFT_result_filtered.csv'
crosslink_file_90_10 = 'C:/Users/jiang/OneDrive/科研/tc/实验记录/2021/20211105 365crosslink/20211115_tangchun_BSA/WXZ_20211108_90_10_HCDFT_result_filtered.csv'
crosslink_file_BSA = 'C:/Users/jiang/OneDrive/科研/tc/实验记录/2021/20211105 365crosslink/20211115_tangchun_BSA/WXZ_20211108_BSA_HCDFT_result_filtered.csv'


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
    pep2_start, pep2_end = find_loc_in_fasta(pep2_fasta, fasta)

    return [[
        pep1_fasta, pep1_num, cross_pos_1, pep1_start, pep1_end,
        pep1_start + int(pep1_num)
    ],
            [
                pep2_fasta, pep2_num, cross_pos_2, pep2_start, pep2_end,
                pep2_start + int(pep2_num)
            ]]


def report_link_pos(ms_file, fasta):
    link_pos_list = []
    with open(ms_file, "r") as f:
        reader = csv.reader(f)
        for target in reader:
            try:
                link_pos_list.append(find_link_pos(target[5], fasta))
            except:
                print("Target link position cannot be analyzed:", target[5])

    return link_pos_list


def cal_repeat_list(list_1, list_2):
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

    return [
        repeat_list, unrepeat_list_1, unrepeat_list_2,
        len(repeat_list) / (len(unrepeat_list_1) + len(repeat_list)),
        len(repeat_list) / (len(unrepeat_list_2) + len(repeat_list))
    ]


def cal_same_link_pos(ms_file_1, ms_file_2, fasta, set=False):
    report_list_1 = report_link_pos(ms_file_1, fasta)
    report_list_2 = report_link_pos(ms_file_2, fasta)
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

    return [
        repeat_list, unrepeat_list_1, unrepeat_list_2, repeat_ratio_1,
        repeat_ratio_2
    ]


def report_animo_ratio(link_list, hard_animo="K"):
    # Input list like [["K","L"],["Y","K"],……]
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

    return cross_dic


def draw_cross_pep_list(dir):
    cross_pep_list = []
    with open(dir) as f:
        reader = csv.reader(f)
        for target in reader:
            if "-" in target[5]:
                cross_pep_list.append(target[5])
    # print(cross_pep_list)
    return cross_pep_list


if __name__ == "__main__":
    repeat_list_3015, unrepeat_list_60_30, unrepeat_list_60_15, ratio_1, ratio_2 = cal_same_link_pos(
        crosslink_file_60_30, crosslink_file_60_15, bsa_fasta, set=True)
    repeat_list_303, unrepeat_list_60_30, unrepeat_list_60_3, ratio_3, ratio_4 = cal_same_link_pos(
        crosslink_file_60_30, crosslink_file_30_10, bsa_fasta, set=True)
    print(repeat_list_3015, repeat_list_303, ratio_1, ratio_2, ratio_3,
          ratio_4)
    print(cal_repeat_list(repeat_list_3015, repeat_list_303))

    cross_list = [[i[0][2], i[1][2]] for i in [
        find_link_pos(j, bsa_fasta)
        for j in draw_cross_pep_list(crosslink_file_60_30)
    ]]
    print(report_animo_ratio(cross_list))

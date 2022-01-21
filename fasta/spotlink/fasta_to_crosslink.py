import csv


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


def report_link_pos(ms_file, fasta, ms_peptide_list=5):
    link_pos_list = []
    with open(ms_file, "r") as f:
        reader = csv.reader(f)
        for target in reader:
            try:
                link_pos_list.append(
                    find_link_pos(target[ms_peptide_list], fasta))
            except:
                print("Target link position cannot be analyzed:", target[5])

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
        return [
            repeat_list, unrepeat_list_1, unrepeat_list_2,
            len(repeat_list) / (len(unrepeat_list_1) + len(repeat_list)),
            len(repeat_list) / (len(unrepeat_list_2) + len(repeat_list))
        ]
    else:
        print(
            len(repeat_list) / (len(unrepeat_list_1) + len(repeat_list)),
            len(repeat_list) / (len(unrepeat_list_2) + len(repeat_list)))
        return [repeat_list]


def fdr2sfdr(file_name):
    return 'site'.join(file_name.split('result'))


def cal_same_link_pos(ms_file_1,
                      ms_file_2,
                      fasta,
                      set=False,
                      simple_mode=False,
                      ms_peptide_list=5,
                      f2sf=False):
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
        print(repeat_ratio_1, repeat_ratio_2)
        return [repeat_list]


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


# Open spotlink file and report ms_peptide
def draw_cross_pep_list(dir):
    cross_pep_list = []
    with open(dir) as f:
        reader = csv.reader(f)
        for target in reader:
            if "-" in target[5]:
                cross_pep_list.append(target[5])
    # print(cross_pep_list)
    return cross_pep_list


# Out dic like {'A': 31, 'C': 26, 'D': 72, 'E': 117, 'F': 11, 'G': 25,……


def report_animo_dic(ms_csv_dir, hard_animo="K"):
    cross_list = [
        [i[0][2], i[1][2]] for i in
        [find_link_pos(j, bsa_fasta) for j in draw_cross_pep_list(ms_csv_dir)]
        if 0 not in [i[0][3], i[0][4], i[1][3], i[1][4]]
    ]
    return report_animo_ratio(cross_list, hard_animo=hard_animo)


if __name__ == "__main__":
    # Def
    bsa_fasta = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA"

    crosslink_file_1108_M10_fDR = 'G:/MSdata/211200BSA/20220120-plink/WXZ_20211108_60_10_HCDFT_result_filtered.csv'
    crosslink_file_1108_M30_fDR = 'G:/MSdata/211200BSA/20220120-plink/WXZ_20211108_60_30_HCDFT_result_filtered.csv'
    crosslink_file_1108_M15_fDR = 'G:/MSdata/211200BSA/20220120-plink/WXZ_20211108_60_15_HCDFT_result_filtered.csv'
    crosslink_file_1108_M3_fDR = 'G:/MSdata/211200BSA/20220120-plink/WXZ_20211108_60_3_HCDFT_result_filtered.csv'
    crosslink_file_1108_N_fDR = 'G:/MSdata/211200BSA/20220120-plink/WXZ_20211108_BSA_HCDFT_result_filtered.csv'

    crosslink_file_0118_M10_fDR = 'G:/MSdata/220118BSA/20220120-plink/JYD_20220114_M10_HCDFT_result_filtered.csv'
    crosslink_file_0118_M30_fDR = 'G:/MSdata/220118BSA/20220120-plink/JYD_20220114_M30_HCDFT_result_filtered.csv'
    crosslink_file_0118_M1_fDR = 'G:/MSdata/220118BSA/20220120-plink/JYD_20220114_M1_HCDFT_result_filtered.csv'
    crosslink_file_0118_S30_fDR = 'G:/MSdata/220118BSA/20220120-plink/JYD_20220114_S30_HCDFT_result_filtered.csv'
    crosslink_file_0118_S15_fDR = 'G:/MSdata/220118BSA/20220120-plink/JYD_20220114_S15_HCDFT_result_filtered.csv'
    crosslink_file_0118_S10_fDR = 'G:/MSdata/220118BSA/20220120-plink/JYD_20220114_S10_HCDFT_result_filtered.csv'
    crosslink_file_0118_S5_fDR = 'G:/MSdata/220118BSA/20220120-plink/JYD_20220114_S5_HCDFT_result_filtered.csv'
    crosslink_file_0118_S1_fDR = 'G:/MSdata/220118BSA/20220120-plink/JYD_20220114_S1_HCDFT_result_filtered.csv'
    crosslink_file_0118_N_fDR = 'G:/MSdata/220118BSA/20220120-plink/JYD_20220114_N_HCDFT_result_filtered.csv'

    repeat_list_s3001 = cal_same_link_pos(crosslink_file_1108_M3_fDR,
                                          crosslink_file_1108_M30_fDR,
                                          bsa_fasta,
                                          set=True,
                                          simple_mode=True,
                                          f2sf=True)[0]

    repeat_list_3010 = cal_same_link_pos(crosslink_file_1108_M10_fDR,
                                         crosslink_file_1108_M30_fDR,
                                         bsa_fasta,
                                         set=True,
                                         simple_mode=True,
                                         f2sf=True)[0]

    repeat_list_N = cal_same_link_pos(crosslink_file_1108_N_fDR,
                                      crosslink_file_1108_M30_fDR,
                                      bsa_fasta,
                                      set=True,
                                      simple_mode=True,
                                      f2sf=True)[0]

    a = cal_repeat_list(repeat_list_3010, repeat_list_s3001,
                        simple_mode=True)[0]
    import pdb_distance_analyze

    for i in a:
        print(
            i,
            pdb_distance_analyze.cal_distance(int(i[0]), int(i[1]),
                                              'G:/MSdata/bsa.pdb'),
            bsa_fasta[int(i[0]) - 1], bsa_fasta[int(i[1]) - 1])

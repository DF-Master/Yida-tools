import sys

sys.path.append('C:/Users/jiang/Documents/GitHub/Yida-tools/fasta/spotlink/')

from fasta_to_crosslink import *
import pdb_distance_analyze
from pymol import cmd

# Def
bsa_fasta = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA"
ub_fasta = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"

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

crosslink_file_0118_UM10_fDR = 'G:/MSdata/220125UBBSA/spotlink/WXZ_20220118_UM10_HCDFT_result_filtered.csv'
crosslink_file_0118_UM1_fDR = 'G:/MSdata/220125UBBSA/spotlink/WXZ_20220118_UM1_HCDFT_result_filtered.csv'
crosslink_file_0118_US30_fDR = 'G:/MSdata/220125UBBSA/spotlink/WXZ_20220118_US30_HCDFT_result_filtered.csv'
crosslink_file_0118_U_fDR = 'G:/MSdata/220125UBBSA/spotlink/WXZ_20220118_U_HCDFT_result_filtered.csv'

crosslink_file_0118_US10_fDR = 'G:/MSdata/220125UBBSA/spotlink/WXZ_20220118_US10_HCDFT_result_filtered.csv'
crosslink_file_0118_US5_fDR = 'G:/MSdata/220125UBBSA/spotlink/WXZ_20220118_US5_HCDFT_result_filtered.csv'
crosslink_file_0118_US1_fDR = 'G:/MSdata/220125UBBSA/spotlink/WXZ_20220118_US1_HCDFT_result_filtered.csv'

# 0118 Time Scale

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
#                                     f2sf=True)[0]

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

#Test 0114 and Test 0118

# repeat_list_01140118_M10 = cal_same_link_pos(crosslink_file_0114_M10_fDR,
#                                              crosslink_file_0118_BM10_fDR,
#                                              bsa_fasta,
#                                              set=True,
#                                              simple_mode=True,
#                                              f2sf=True)[0]

# repeat_list_01140118_S30 = cal_same_link_pos(crosslink_file_0114_S30_fDR,
#                                              crosslink_file_0118_BS30_fDR,
#                                              bsa_fasta,
#                                              set=True,
#                                              simple_mode=True,
#                                              f2sf=True)[0]

# repeat_list_01140118_S10 = cal_same_link_pos(crosslink_file_0114_S10_fDR,
#                                              crosslink_file_0118_BS10_fDR,
#                                              bsa_fasta,
#                                              set=True,
#                                              simple_mode=True,
#                                              f2sf=True)[0]

# repeat_list_01140118_S5 = cal_same_link_pos(crosslink_file_0114_S5_fDR,
#                                             crosslink_file_0118_BS5_fDR,
#                                             bsa_fasta,
#                                             set=True,
#                                             simple_mode=True,
#                                             f2sf=True)[0]

# repeat_list_01140118_N = cal_same_link_pos(crosslink_file_0114_N_fDR,
#                                            crosslink_file_0118_B_fDR,
#                                            bsa_fasta,
#                                            set=True,
#                                            simple_mode=True,
#                                            f2sf=True)[0]

# a = cal_repeat_list(repeat_list_01140118_M10,
#                     repeat_list_01140118_S30,
#                     simple_mode=True)[0]
# b = cal_repeat_list(repeat_list_01140118_S10,
#                     repeat_list_01140118_S5,
#                     simple_mode=True)[0]

# c = cal_repeat_list(a, b, simple_mode=True)[0]

# report_animo_ratio(c, pos_form=True, fasta_file=bsa_fasta, tableform=True)

# num = 0
# for i in a:
#     dis = pdb_distance_analyze.cal_distance(int(i[0]), int(i[1]),
#                                             'G:/MSdata/bsa.pdb')
#     print(i, dis, bsa_fasta[int(i[0]) - 1], bsa_fasta[int(i[1]) - 1])
#     # cmd.distance(str(i[0]) + '-' + str(i[1]),
#     #              "(/bsa//A/" + str(i[0]) + "/CA)",
#     #              "( /bsa//A/" + str(i[1]) + "/CA)",
#     #              mode=0)
#     if dis >= 15:
#         num += 1
# print('Over 15 A:', num)

# Test 0118 Ub

# repeat_list_0118_UM1 = cal_same_link_pos(crosslink_file_0118_UM10_fDR,
#                                          crosslink_file_0118_UM1_fDR,
#                                          ub_fasta,
#                                          set=True,
#                                          simple_mode=True,
#                                          f2sf=True)[0]
# repeat_list_0118_US30 = cal_same_link_pos(crosslink_file_0118_UM10_fDR,
#                                           crosslink_file_0118_US30_fDR,
#                                           ub_fasta,
#                                           set=True,
#                                           simple_mode=True,
#                                           f2sf=True)[0]
# # repeat_list_0118_U = cal_same_link_pos(crosslink_file_0118_UM10_fDR,
# #                                        crosslink_file_0118_U_fDR,
# #                                        ub_fasta,
# #                                        set=True,
# #                                        simple_mode=True,
# #                                        f2sf=True)[0]
# repeat_list_0118_US10 = cal_same_link_pos(crosslink_file_0118_UM10_fDR,
#                                           crosslink_file_0118_US10_fDR,
#                                           ub_fasta,
#                                           set=True,
#                                           simple_mode=True,
#                                           f2sf=True)[0]
# repeat_list_0118_US5 = cal_same_link_pos(crosslink_file_0118_UM10_fDR,
#                                          crosslink_file_0118_US5_fDR,
#                                          ub_fasta,
#                                          set=True,
#                                          simple_mode=True,
#                                          f2sf=True)[0]
# repeat_list_0118_US1 = cal_same_link_pos(crosslink_file_0118_UM10_fDR,
#                                          crosslink_file_0118_US1_fDR,
#                                          ub_fasta,
#                                          set=True,
#                                          simple_mode=True,
#                                          f2sf=True)[0]

# US30 = cal_repeat_list(repeat_list_0118_US30,
#                        repeat_list_0118_UM1,
#                        simple_mode=True)[0]
# US10 = cal_repeat_list(repeat_list_0118_US10,
#                        repeat_list_0118_UM1,
#                        simple_mode=True)[0]
# US5 = cal_repeat_list(repeat_list_0118_US5,
#                       repeat_list_0118_UM1,
#                       simple_mode=True)[0]
# US1 = cal_repeat_list(repeat_list_0118_US1,
#                       repeat_list_0118_UM1,
#                       simple_mode=True)[0]
# allin = US30
# for i in [US10, US5, US1]:
#     allin = cal_repeat_list(allin, i, simple_mode=True)[0]

# report_animo_ratio(allin, pos_form=True, fasta_file=ub_fasta, tableform=True)

# num = 0
# for i in allin:
#     dis = pdb_distance_analyze.cal_distance(
#         int(i[0]), int(i[1]), 'G:/MSdata/220125UBBSA/mono_ub.pdb')
#     print(i, dis, ub_fasta[int(i[0]) - 1], ub_fasta[int(i[1]) - 1])
#     # cmd.distance(str(i[0]) + '-' + str(i[1]),
#     #              "(/mono_ub//A/" + str(i[0]) + "/CA)",
#     #              "( /mono_ub//A/" + str(i[1]) + "/CA)",
#     #              mode=0)
#     if dis >= 15:
#         num += 1
# print('Over 15 A:', num)

# Selection

repeat_list_01140118_M10 = cal_same_link_pos(crosslink_file_0114_M10_fDR,
                                             crosslink_file_0118_BM10_fDR,
                                             bsa_fasta,
                                             set=True,
                                             simple_mode=True,
                                             f2sf=True)[0]

repeat_list_01140118_S30 = cal_same_link_pos(crosslink_file_0114_S30_fDR,
                                             crosslink_file_0118_BS30_fDR,
                                             bsa_fasta,
                                             set=True,
                                             simple_mode=True,
                                             f2sf=True)[0]

list_0118_M10 = cal_same_link_pos(crosslink_file_0118_BM10_fDR,
                                  crosslink_file_0118_BM10_fDR,
                                  bsa_fasta,
                                  set=True,
                                  simple_mode=True,
                                  f2sf=True)[0]

list_0118_S5 = cal_same_link_pos(crosslink_file_0118_BS5_fDR,
                                 crosslink_file_0118_BS5_fDR,
                                 bsa_fasta,
                                 set=True,
                                 simple_mode=True,
                                 f2sf=True)[0]
num = 0
for i in repeat_list_01140118_M10:
    dis = pdb_distance_analyze.cal_distance(int(i[0]), int(i[1]),
                                            'G:/MSdata/bsa.pdb')
    print(i, dis, bsa_fasta[int(i[0]) - 1], bsa_fasta[int(i[1]) - 1])
    if dis >= 15:
        num += 1
print('Over 15 A:', num)